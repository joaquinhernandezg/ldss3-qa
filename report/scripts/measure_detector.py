"""
Detector measurements for the LDSS3-C commissioning report.

Works directly on the raw amplifier files (``*c1.fits``/``*c2.fits``), using the
PypeIt LDSS3 reader only to assemble the two halves of each exposure, so that
the numbers are independent of any reduction choice.

Measured quantities, per observing epoch and per amplifier:

* bias (overscan) level and its row-to-row stability;
* read noise, from the robust scatter of the difference of consecutive bias
  frames, and independently from the overscan of each bias;
* gain, from a photon-transfer analysis of pairs of flat-field exposures
  (variance of the pair difference against the mean signal);
* dark current, from long dark exposures bracketed by biases;
* the column profile of a combined bias, used to define the bad-column mask;
* the header keywords ``SPEED``, ``EGAIN`` and ``ENOISE``.

Usage::

    python measure_detector.py <raw_root> <results_dir>

``<raw_root>`` holds one directory per night of raw data.  Results are written
as ``detector.json`` plus ``bias_profiles.npz`` and ``ptc.npz`` for plotting.
"""
import sys
import json
from pathlib import Path
from collections import defaultdict

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from pypeit.spectrographs.util import load_spectrograph

SPEC = load_spectrograph('magellan_ldss3')

# Layout of the assembled raw frame returned by get_rawimage: (spectral,
# spatial) = (4096, 2304); amplifier 1 overscan in columns 0-127, data in
# 128-1151; amplifier 2 data in 1152-2175, overscan in 2176-2303.
AMP_DATA = {1: slice(128, 1152), 2: slice(1152, 2176)}
AMP_OSCAN = {1: slice(0, 128), 2: slice(2176, 2304)}
# Region used for statistics: avoid the first/last 100 rows, the outer 20
# columns of each amplifier and the columns next to the seam.
ROWS = slice(100, 3996)
COLS = {1: slice(20, 1004), 2: slice(20, 1004)}

# Bad-column ranges currently in the LDSS3 bad-pixel mask (trimmed frame)
BPM_RANGES = [(0, 11), (443, 443), (608, 608), (1413, 1413), (1492, 1494),
              (1549, 1551), (1602, 1606), (1635, 1640), (1688, 1690),
              (1695, 1695), (1999, 2000), (2034, 2047)]


def robust_sigma(x):
    """Standard deviation after iterative 4-sigma clipping.

    The raw data are integers, so median/MAD-based estimators are quantised
    (to half-integer ADU) and are not used for the noise measurements."""
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    _, _, std = sigma_clipped_stats(x, sigma=4, maxiters=10, cenfunc='median', stdfunc='std')
    return float(std)


def read(path):
    """Assembled, overscan-subtracted frame in ADU, trimmed to 4096 x 2048,
    plus the per-amplifier overscan arrays and header."""
    _, arr, hdu, exptime, _, _ = SPEC.get_rawimage(str(path), 1)
    hdr = hdu[0].header
    arr = arr.astype(float)
    trimmed = np.empty((arr.shape[0], 2048))
    oscan = {}
    for amp in (1, 2):
        # Skip the first few overscan columns, which carry charge from the
        # last data column.
        osc = arr[:, AMP_OSCAN[amp]][:, 10:-5]
        oscan[amp] = osc
        level = np.median(osc, axis=1)
        data = arr[:, AMP_DATA[amp]] - level[:, None]
        trimmed[:, (amp - 1) * 1024:amp * 1024] = data
    return trimmed, oscan, hdr, exptime


def amp_view(img, amp):
    return img[:, (amp - 1) * 1024:amp * 1024][ROWS, COLS[amp]]


def bpm_columns():
    good = np.ones(2048, dtype=bool)
    for c1, c2 in BPM_RANGES:
        good[c1:c2 + 1] = False
    return good


def inventory(raw_root):
    """Header survey of every amplifier-1 file."""
    rows = []
    for night in sorted(Path(raw_root).iterdir()):
        for f in sorted(night.glob('ccd*c1.fits')):
            h = fits.getheader(f)
            rows.append(dict(file=str(f), night=night.name, object=str(h.get('OBJECT', '')),
                             exptime=float(h.get('EXPTIME', -1)), grism=h.get('GRISM'),
                             filter=h.get('FILTER'), aperture=h.get('APERTURE'),
                             speed=h.get('SPEED'), egain=h.get('EGAIN'),
                             enoise=h.get('ENOISE'), date=h.get('UT-DATE'),
                             time=h.get('UT-TIME'),
                             egain2=fits.getheader(str(f).replace('c1.fits', 'c2.fits')).get('EGAIN'),
                             enoise2=fits.getheader(str(f).replace('c1.fits', 'c2.fits')).get('ENOISE')))
    return rows


def is_bias(r):
    o = r['object'].lower()
    return r['exptime'] == 0 and ('bias' in o or 'zero' in o)


def is_dark(r):
    return 'dark' in r['object'].lower() and r['exptime'] > 0


def epoch_of(r):
    return r['date'][:7]          # YYYY-MM


def measure_bias(files):
    """Read noise and bias statistics from a set of bias frames."""
    frames, oscans, hdrs = [], [], []
    for f in files:
        img, osc, hdr, _ = read(f)
        frames.append(img)
        oscans.append(osc)
        hdrs.append(hdr)
    gain = {1: float(hdrs[0]['EGAIN']),
            2: float(fits.getheader(files[0].replace('c1.fits', 'c2.fits'))['EGAIN'])}
    out = {}
    for amp in (1, 2):
        # Frame differences: consecutive pairs
        diff_rn = []
        for a, b in zip(frames[:-1], frames[1:]):
            d = amp_view(a, amp) - amp_view(b, amp)
            diff_rn.append(robust_sigma(d) / np.sqrt(2))
        diff_rn = np.array(diff_rn)
        # Overscan scatter, after removing the row median
        osc_rn = np.array([robust_sigma((o[amp] - np.median(o[amp], axis=1)[:, None])[ROWS])
                           for o in oscans])
        # Raw overscan level: from the un-subtracted data
        out[amp] = dict(
            n_frames=len(files),
            gain_hdr=gain[amp],
            ron_adu=float(np.mean(diff_rn)), ron_adu_std=float(np.std(diff_rn)),
            ron_e=float(np.mean(diff_rn) * gain[amp]),
            ron_e_std=float(np.std(diff_rn) * gain[amp]),
            ron_oscan_e=float(np.mean(osc_rn) * gain[amp]),
            # residual 2D structure left in the overscan-subtracted bias
            bias_resid_adu=float(np.median(np.median(np.array([amp_view(f, amp) for f in frames]), axis=0))),
        )
    # Column profile: clipped mean over frames and rows.  Medians of the
    # integer-valued raw data would be quantised to half-integer ADU.
    stack = np.array(frames)[:, ROWS, :]
    profile = sigma_clipped_stats(stack.reshape(-1, stack.shape[-1]), sigma=3, maxiters=5,
                                  axis=0)[0]
    return out, profile, [dict(speed=h.get('SPEED'), enoise=h.get('ENOISE')) for h in hdrs[:1]]


def overscan_levels(files):
    lv = defaultdict(list)
    for f in files:
        _, arr, _, _, _, _ = SPEC.get_rawimage(str(f), 1)
        for amp in (1, 2):
            lv[amp].append(float(np.median(arr[:, AMP_OSCAN[amp]][:, 10:-5])))
    return {a: (float(np.mean(v)), float(np.std(v))) for a, v in lv.items()}


def ptc(pairs, ron_adu):
    """Photon-transfer gain from pairs of flats.

    For each pair the second frame is scaled to the median of the first
    (lamp drift), the difference is formed, and pixels are binned by mean
    signal.  In each bin var(D)/2 - RN^2 = S/g, so g follows from a linear fit.
    """
    res = {}
    curves = {}
    for amp in (1, 2):
        S_all, V_all = [], []
        for fa, fb in pairs:
            a = amp_view(read(fa)[0], amp)
            b = amp_view(read(fb)[0], amp)
            lit = a > 200
            scale = np.median(a[lit]) / np.median(b[lit])
            b = b * scale
            s = 0.5 * (a + b)
            d = a - b
            # Keep only pixels where the illumination is locally smooth, so that
            # small shifts of the slit image or of the lamp pattern between the
            # two exposures do not add variance (slit edges, bridges).
            gx = np.zeros_like(s)
            gy = np.zeros_like(s)
            gx[:, 1:-1] = np.abs(s[:, 2:] - s[:, :-2]) / 2
            gy[1:-1, :] = np.abs(s[2:, :] - s[:-2, :]) / 2
            lit &= (gx < 0.005 * s) & (gy < 0.005 * s)
            # bins of signal; within a bin use a robust variance
            edges = np.quantile(s[lit], np.linspace(0.02, 0.98, 30))
            idx = np.digitize(s, edges)
            for i in range(1, len(edges)):
                m = (idx == i) & lit & (s < 50000)
                if m.sum() < 2000:
                    continue
                S_all.append(np.median(s[m]))
                V_all.append(robust_sigma(d[m]) ** 2 / 2 - ron_adu[amp] ** 2)
        S_all, V_all = np.array(S_all), np.array(V_all)
        good = (S_all > 500) & (S_all < 40000)
        # Fit V = S/g through the origin, with 3-sigma clipping
        for _ in range(3):
            inv_g = np.sum(S_all[good] * V_all[good]) / np.sum(S_all[good] ** 2)
            r = V_all - inv_g * S_all
            good &= np.abs(r) < 3 * robust_sigma(r[good])
        # Bootstrap for the uncertainty
        rng = np.random.default_rng(1234)
        boots = []
        gi = np.where(good)[0]
        for _ in range(500):
            j = rng.choice(gi, gi.size)
            boots.append(1 / (np.sum(S_all[j] * V_all[j]) / np.sum(S_all[j] ** 2)))
        res[amp] = dict(gain=float(1 / inv_g), gain_err=float(np.std(boots)),
                        n_points=int(good.sum()), n_pairs=len(pairs))
        curves[amp] = np.vstack([S_all, V_all, good])
    return res, curves


def dark_current(dark_files, bias_combined, gain):
    out = {}
    for amp in (1, 2):
        rates = []
        for f in dark_files:
            img, _, _, t = read(f)
            d = amp_view(img - bias_combined, amp)
            mean, _, _ = sigma_clipped_stats(d, sigma=3, maxiters=10)
            rates.append(mean * gain[amp] / t * 3600.)
        out[amp] = dict(rate_e_hr=float(np.mean(rates)), std=float(np.std(rates)),
                        n=len(dark_files),
                        exptimes=sorted(set(float(fits.getheader(f)['EXPTIME']) for f in dark_files)))
    return out


def main(raw_root, results_dir):
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    inv = inventory(raw_root)

    # --- header survey --------------------------------------------------------
    survey = defaultdict(set)
    for r in inv:
        survey[epoch_of(r)].add((r['speed'], r['egain'], r['enoise'], r['egain2'], r['enoise2']))
    header_survey = {k: sorted([list(v) for v in s]) for k, s in survey.items()}

    # --- biases by epoch ------------------------------------------------------
    biases = defaultdict(list)
    for r in inv:
        if is_bias(r):
            biases[(epoch_of(r), r['night'])].append(r['file'])
    bias_res, profiles, levels = {}, {}, {}
    for (ep, night), files in sorted(biases.items()):
        if len(files) < 3:
            continue
        key = f'{ep} ({night})'
        print(f'bias {key}: {len(files)} frames')
        out, prof, _ = measure_bias(sorted(files))
        bias_res[key] = out
        profiles[key] = prof
        levels[key] = overscan_levels(sorted(files))

    # --- gain from flat pairs ---------------------------------------------------
    # Pairs of consecutive flats with identical set-up and exposure time
    flats = [r for r in inv if any(k in r['object'].lower() for k in ('flat', 'qh', ' ff'))
             and r['exptime'] > 0]
    groups = defaultdict(list)
    for r in flats:
        groups[(r['night'], r['grism'], r['aperture'], r['filter'], r['exptime'])].append(r['file'])
    ptc_res, ptc_curves = {}, {}
    for key, files in sorted(groups.items()):
        files = sorted(files)
        if len(files) < 2:
            continue
        pairs = list(zip(files[0::2], files[1::2]))[:4]
        # read noise in ADU for this epoch: use the nearest bias measurement,
        # falling back on the published Fast values
        ron_adu = {1: 7.0 / 1.65, 2: 7.2 / 1.47}
        name = f'{key[0]} {key[1]} {key[2]} {key[3]} {key[4]:g}s'
        print(f'ptc {name}: {len(pairs)} pairs')
        r, c = ptc(pairs, ron_adu)
        ptc_res[name] = r
        ptc_curves[name] = c

    # --- dark current -----------------------------------------------------------
    darks = defaultdict(list)
    for r in inv:
        if is_dark(r):
            darks[r['night']].append(r['file'])
    dark_res = {}
    for night, files in darks.items():
        bfiles = [f for (ep, n), fl in biases.items() if n == night for f in fl]
        if len(bfiles) < 3:
            continue
        bias_comb = np.median(np.array([read(f)[0] for f in sorted(bfiles)]), axis=0)
        h1 = fits.getheader(files[0])
        h2 = fits.getheader(files[0].replace('c1.fits', 'c2.fits'))
        dark_res[night] = dark_current(sorted(files), bias_comb,
                                       {1: float(h1['EGAIN']), 2: float(h2['EGAIN'])})
        print(f'dark {night}: {dark_res[night]}')

    # --- bad columns from the bias profile ----------------------------------------
    good = bpm_columns()
    flagged = {}
    for key, prof in profiles.items():
        resid = prof - np.array([np.median(prof[max(0, i - 25):i + 25]) for i in range(prof.size)])
        sig = robust_sigma(resid[good])
        flagged[key] = [int(i) for i in np.where(np.abs(resid) > 8 * sig)[0]]

    out = dict(header_survey=header_survey, bias=bias_res, overscan_level=levels,
               ptc=ptc_res, dark=dark_res, bad_column_candidates=flagged,
               bpm_ranges=BPM_RANGES)
    with open(results_dir / 'detector.json', 'w') as fh:
        json.dump(out, fh, indent=1, default=str)
    np.savez(results_dir / 'bias_profiles.npz', **{k: v for k, v in profiles.items()})
    np.savez(results_dir / 'ptc.npz', **{f'{k}|amp{a}': c[a] for k, c in ptc_curves.items()
                                          for a in (1, 2)})
    print(json.dumps(out, indent=1, default=str)[:4000])


if __name__ == '__main__':
    main(*sys.argv[1:3])
