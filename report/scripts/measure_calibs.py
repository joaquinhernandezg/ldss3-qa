"""
Calibration and reduction diagnostics for the commissioning report.

For each reduced dataset (see ``style.DATASETS``) this collects, from the
PypeIt output files only:

* slit tracing: number of slits, edges and widths at mid-detector, mask flags;
* wavelength calibration: per-slit number of lines, fit RMS, dispersion,
  coverage, arc-line FWHM and resolving power, and the per-line residuals;
* spectral tilts: RMS of the 2D tilt fit to the traced arc lines;
* flat field: the flux step at the amplifier seam before and after
  normalisation, and the peak flat level per slit;
* science frames: sky-line flexure shifts, the continuity of the wavelength
  image across the seam, the noise-model check sigma(chi), and the S/N of
  extracted objects.

Results go to ``results/calibs.json``.
"""
import sys
import json
import warnings
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from ldss3_qa import step_across                                  # noqa: E402
from style import REDUX, RESULTS, DATASETS                        # noqa: E402

from pypeit.slittrace import SlitTraceSet                         # noqa: E402
from pypeit.wavecalib import WaveCalib                            # noqa: E402
from pypeit.wavetilts import WaveTilts                            # noqa: E402
from pypeit.flatfield import FlatImages                           # noqa: E402
from pypeit.spec2dobj import AllSpec2DObj                         # noqa: E402
from pypeit.specobjs import SpecObjs                              # noqa: E402

SEAM = 1024
PLATESCALE = 0.189


def robust_sigma(x):
    """Standard deviation after iterative 5-sigma clipping."""
    from astropy.stats import sigma_clipped_stats
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return float(sigma_clipped_stats(x, sigma=5, maxiters=10)[2])


def calib_files(d, prefix):
    return sorted((d / 'Calibrations').glob(f'{prefix}_*'))


def slit_info(path):
    st = SlitTraceSet.from_file(str(path), chk_version=False)
    mid = st.nspec // 2
    left, right, _ = st.select_edges(initial=True)
    rows = []
    for i in range(st.nslits):
        flags = st.bitmask.flagged_bits(st.mask[i]) if st.mask[i] else []
        rows.append(dict(spat_id=int(st.spat_id[i]), left=float(left[mid, i]),
                         right=float(right[mid, i]),
                         width_pix=float(right[mid, i] - left[mid, i]),
                         width_arcsec=float((right[mid, i] - left[mid, i]) * PLATESCALE),
                         crosses_seam=bool(left[mid, i] < SEAM < right[mid, i]),
                         flags=list(flags)))
    return dict(nslits=int(st.nslits), slits=rows)


def wave_info(path):
    wc = WaveCalib.from_file(str(path), chk_version=False)
    rows, lines = [], []
    for i, f in enumerate(wc.wv_fits):
        if f is None or f.pypeitfit is None or f.wave_soln is None or f.wave_fit is None:
            rows.append(dict(spat_id=int(wc.spat_ids[i]), ok=False))
            continue
        ws = f.wave_soln
        gpm = f.pypeitfit.gpm.astype(bool)
        model = f.pypeitfit.eval(f.pixel_fit / f.xnorm)
        resid = f.wave_fit - model
        disp = float(np.abs(np.median(np.diff(ws))))
        ions = f.ions if f.ion_bits is not None else ['?'] * f.wave_fit.size
        rows.append(dict(
            spat_id=int(wc.spat_ids[i]), ok=True, nlin=int(gpm.sum()),
            nlin_total=int(f.wave_fit.size), rms_pix=float(f.rms),
            rms_ang=float(np.sqrt(np.mean(resid[gpm] ** 2))),
            disp=disp, cen_wave=float(f.cen_wave),
            wmin=float(ws.min()), wmax=float(ws.max()),
            wline_min=float(f.wave_fit[gpm].min()), wline_max=float(f.wave_fit[gpm].max()),
            fwhm_pix=float(f.fwhm), fwhm_ang=float(f.fwhm * disp),
            R=float(f.cen_wave / (f.fwhm * disp)),
            order=int(f.pypeitfit.order[0])))
        for w, p, r, g, ion in zip(f.wave_fit, f.pixel_fit, resid, gpm, ions):
            lines.append(dict(spat_id=int(wc.spat_ids[i]), wave=float(w), pix=float(p),
                              resid=float(r), used=bool(g), ion=str(ion)))
    return dict(lamps=str(wc.lamps), slits=rows), lines


def tilt_info(path):
    wt = WaveTilts.from_file(str(path), chk_version=False)
    out = dict(spat_order=[int(x) for x in np.atleast_1d(wt.spat_order)],
               spec_order=[int(x) for x in np.atleast_1d(wt.spec_order)])
    tt = wt.tilt_traces
    if tt is not None:
        out['tilt_trace_columns'] = list(tt.colnames) if hasattr(tt, 'colnames') else None
    return out


def flat_info(path, slits_path):
    fi = FlatImages.from_file(str(path), chk_version=False)
    out = {}
    raw = fi.pixelflat_raw
    st = SlitTraceSet.from_file(str(slits_path), chk_version=False)
    if raw is not None:
        # The step is measured inside a slit that spans the seam, using at most
        # the part of the slit on each side (minus a few pixels at the edges).
        mid = raw.shape[0] // 2
        room = [min(SEAM - st.left_init[mid, i], st.right_init[mid, i] - SEAM)
                for i in range(st.nslits)]
        i = int(np.argmax(room))
        win = int(room[i]) - 8 - 3
        out['seam_slit'] = int(st.spat_id[i]) if win >= 12 else None
        if win >= 12:
            rows = slice(mid - 1000, mid + 1000)
            j, _, _ = step_across(np.median(raw[rows], axis=0), win=min(win, 110))
            out['seam_step_raw_pct'] = None if j is None else float(j)
            if fi.pixelflat_norm is not None:
                jn, _, _ = step_across(np.median(fi.pixelflat_norm[rows], axis=0),
                                       win=min(win, 110))
                out['seam_step_norm_pct'] = None if jn is None else float(jn)
        peaks = []
        for i in range(st.nslits):
            l, r = int(max(st.left_init[mid, i], 0)), int(min(st.right_init[mid, i], raw.shape[1]))
            sub = raw[mid - 400:mid + 400, l + 2:r - 2]
            peaks.append(float(np.nanpercentile(sub, 99)) if sub.size else np.nan)
        out['flat_p99_e'] = peaks
        if fi.pixelflat_norm is not None:
            norm = fi.pixelflat_norm
            good = (norm > 0.5) & (norm < 1.5) & (norm != 1.0)
            out['pixelflat_rms'] = robust_sigma(norm[good])
    return out


def seam_wave_jump(s2):
    """Wavelength discontinuity across the seam, in Angstrom.  For every row,
    extrapolate the wavelength image linearly from each side to the seam."""
    w = s2.waveimg
    jumps = []
    for row in range(200, w.shape[0] - 200, 25):
        wl, wr = w[row, SEAM - 12:SEAM - 2], w[row, SEAM + 2:SEAM + 12]
        if np.any(wl <= 0) or np.any(wr <= 0):
            continue
        xl, xr = np.arange(SEAM - 12, SEAM - 2), np.arange(SEAM + 2, SEAM + 12)
        el = np.polyval(np.polyfit(xl, wl, 1), SEAM - 0.5)
        er = np.polyval(np.polyfit(xr, wr, 1), SEAM - 0.5)
        jumps.append(er - el)
    if not jumps:
        return None
    return dict(median=float(np.median(jumps)), max_abs=float(np.max(np.abs(jumps))),
                n_rows=len(jumps))


ENOISE_CARD = (4.67, 5.06)       # value written by the instrument in every frame
HIST_EDGES = np.linspace(-8, 8, 161)
HISTS = {'published': [], 'enoise': []}


def noise_chi(s2, exclude=None, keep_hist=True):
    """sigma of chi = (data - sky - object) * sqrt(ivar) over good pixels.

    The reduction uses the published read noise for the readout mode.  The
    effect of using the ENOISE header card instead is evaluated by swapping the
    read-noise term in the variance, amplifier by amplifier."""
    resid = s2.sciimg - s2.skymodel - s2.objmodel
    gd = (s2.bpmmask.mask == 0) & np.isfinite(resid) & (s2.waveimg > 0) & (s2.ivarmodel > 0)
    var = np.zeros_like(resid)
    var[gd] = 1.0 / s2.ivarmodel[gd]
    # Sky-dominated pixels only: residuals of bright objects (standards,
    # alignment stars) are dominated by profile-model errors, not by noise.
    gd &= np.abs(s2.objmodel) < 0.1 * np.sqrt(np.where(gd, var, 1.0))
    if exclude is not None:
        gd &= ~exclude
    rn = np.atleast_1d(s2.detector['ronoise']).astype(float)
    amp = np.ones(resid.shape, dtype=int)
    amp[:, SEAM:] = 2
    var_en = var.copy()
    for a in (1, 2):
        m = gd & (amp == a)
        var_en[m] = var[m] - rn[a - 1] ** 2 + ENOISE_CARD[a - 1] ** 2
    out = dict(npix=int(gd.sum()), ronoise=[float(r) for r in rn],
               gain=[float(g) for g in np.atleast_1d(s2.detector['gain'])])
    # Sky-subtraction residuals at bright sky lines: chi over the 10 per cent of
    # sky pixels with the brightest sky model, against the remaining 90.
    sk = s2.skymodel[gd]
    thr = np.percentile(sk, 90)
    chi_pub = resid[gd] / np.sqrt(var[gd])
    out['sigma_brightsky'] = robust_sigma(chi_pub[sk > thr])
    out['sigma_faintsky'] = robust_sigma(chi_pub[sk <= thr])
    for key, v in (('published', var), ('enoise', var_en)):
        ok = gd & (v > 0)
        chi = resid[ok] / np.sqrt(v[ok])
        out[f'sigma_{key}'] = robust_sigma(chi)
        if keep_hist:
            HISTS[key].append(np.histogram(chi, bins=HIST_EDGES)[0])
    return out


OI_VAC = [5578.887, 6302.046, 6365.536]   # [OI] 5577.339, 6300.304, 6363.776 A (air)


def _gauss(x, a, m, sig, b):
    return a * np.exp(-0.5 * ((x - m) / sig) ** 2) + b


def sky_line_offsets(s2, slitmask_ids):
    """Offset (measured - reference, in A) of sky emission lines in the final,
    flexure-corrected wavelength frame, per slit.

    Each line is fitted with a Gaussian plus constant.  Fits whose width departs
    from the median width by more than 40 per cent (blends, or no line) are
    discarded.  References: the [OI] night-sky lines and the brightest lines of
    OH_LDSS3_vac that have no other listed line within 8 A."""
    from scipy.optimize import curve_fit
    from pypeit.core.wavecal import waveio
    oh = waveio.load_line_lists(['OH_LDSS3_vac'])[0]
    ohw, oha = np.asarray(oh['wave']), np.asarray(oh['amplitude'])
    iso = np.array([np.min(np.abs(np.delete(ohw, i) - w)) > 8 for i, w in enumerate(ohw)])
    ref = np.concatenate([ohw[iso & (oha > np.percentile(oha, 50))], OI_VAC])
    vc = s2.vel_corr if s2.vel_corr not in (None, 0) else 1.0
    res = {}
    for sid in np.unique(slitmask_ids):
        if sid <= 0:
            continue
        cols = np.where((slitmask_ids == sid).any(axis=0))[0]
        c = int(np.median(cols))
        # The saved wavelength image includes the heliocentric correction; sky
        # lines are at rest in the topocentric frame.
        w = s2.waveimg[:, c] / vc
        sky = np.median(s2.skymodel[:, max(cols[0] + 3, c - 3):min(cols[-1] - 3, c + 4)], axis=1)
        good = w > 0
        if good.sum() < 1000:
            continue
        w, sky = w[good], sky[good]
        fits = []
        for lam in ref:
            k = int(np.argmin(np.abs(w - lam)))
            if k < 10 or k > w.size - 10 or abs(w[k] - lam) > 2:
                continue
            sl = slice(k - 6, k + 7)
            y = sky[sl]
            try:
                p, _ = curve_fit(_gauss, w[sl], y, p0=[y.max() - y.min(), lam, 2.5, y.min()],
                                 maxfev=2000)
            except Exception:                                       # noqa: BLE001
                continue
            if p[0] <= 0 or abs(p[1] - lam) > 4:
                continue
            fits.append((lam, p[1] - lam, abs(p[2])))
        if len(fits) < 2:
            continue
        f = np.array(fits)
        msig = np.median(f[:, 2])
        ok = np.abs(f[:, 2] / msig - 1) < 0.4
        d = f[ok, 1]
        if d.size < 2:
            continue
        ok2 = np.abs(d - np.median(d)) < max(3 * robust_sigma(d), 0.3) if d.size > 3 \
            else np.ones(d.size, bool)
        d = d[ok2]
        res[int(sid)] = dict(n=int(d.size), median=float(np.median(d)),
                             scatter=float(np.std(d)), err=float(np.std(d) / np.sqrt(d.size)),
                             lines=[float(x) for x in f[ok][ok2][:, 0]],
                             offsets=[float(x) for x in f[ok][ok2][:, 1]],
                             disp=float(np.median(np.diff(w))))
    return res


def box_slits(ds):
    """Alignment boxes: slits whose arc lines are much broader than the
    median of the mask (the boxes are wide in the dispersion direction) or
    whose wavelength calibration failed."""
    wp = sorted((REDUX / ds / 'Calibrations').glob('WaveCalib_*.fits'))
    if not wp:
        return set()
    wc = WaveCalib.from_file(str(wp[0]), chk_version=False)
    fw = np.array([f.fwhm if (f is not None and f.pypeitfit is not None) else np.nan
                   for f in wc.wv_fits])
    # Reference: the narrow slits.  A mask can have as many holes as slits, so
    # the median is not a safe reference.
    ref = np.nanpercentile(fw, 20)
    return {int(s) for s, f in zip(wc.spat_ids, fw) if not np.isfinite(f) or f > 1.5 * ref}


def science_info(d):
    """Diagnostics for every reduced exposure.  Standards are reduced without a
    global sky model, so only frames typed 'science' enter the noise test."""
    from pypeit.inputfiles import PypeItFile
    tbl = PypeItFile.from_file(str(d / f'{d.name}.pypeit')).data
    ftype = {str(f).replace('.fits', ''): str(t) for f, t in zip(tbl['filename'], tbl['frametype'])}
    out = []
    boxes = box_slits(d.name)
    for f2 in sorted((d / 'Science').glob('spec2d_*.fits')):
        a = AllSpec2DObj.from_fits(str(f2), chk_version=False)
        s2 = a[a.detectors[0]]
        slitmask = s2.slits.slit_img(initial=True, flexure=s2.sci_spat_flexure) \
            if hasattr(s2.slits, 'slit_img') else None
        exclude = np.isin(slitmask, list(boxes)) if slitmask is not None else None
        stem = f2.name.split('_')[1].split('-')[0]
        is_sci = 'science' in ftype.get(stem, '')
        rec = dict(file=f2.name, frametype=ftype.get(stem, ''), seam_wave_jump=seam_wave_jump(s2),
                   noise=noise_chi(s2, exclude=exclude, keep_hist=is_sci), boxes=sorted(boxes),
                   sky_lines=sky_line_offsets(s2, slitmask) if slitmask is not None else None)
        f1 = f2.with_name(f2.name.replace('spec2d_', 'spec1d_'))
        if f1.exists():
            so = SpecObjs.from_fitsfile(str(f1), chk_version=False)
            objs = []
            for s in so:
                flex = getattr(s, 'FLEX_SHIFT_TOTAL', None)
                sn = None
                if s.OPT_COUNTS is not None:
                    m = s.OPT_MASK & (s.OPT_COUNTS_IVAR > 0)
                    sn = float(np.median(s.OPT_COUNTS[m] * np.sqrt(s.OPT_COUNTS_IVAR[m]))) \
                        if m.any() else None
                objs.append(dict(name=s.NAME, slit=int(s.SLITID), spat_pix=float(s.SPAT_PIXPOS),
                                 in_box=int(s.SLITID) in boxes,
                                 flex_shift=None if flex is None else float(flex),
                                 sn_median=sn, fwhm_arcsec=float(s.FWHM) if s.FWHM is not None else None))
            rec['objects'] = objs
        out.append(rec)
    return out


def main():
    RESULTS.mkdir(exist_ok=True)
    res, all_lines = {}, {}
    for name in DATASETS:
        d = REDUX / name
        if not d.is_dir():
            continue
        print(f'== {name}')
        r = {}
        for sp in calib_files(d, 'Slits'):
            key = sp.name.replace('Slits_', '').replace('.fits.gz', '')
            c = dict(slits=slit_info(sp))
            wp = d / 'Calibrations' / f'WaveCalib_{key}.fits'
            if wp.exists():
                c['wave'], all_lines[f'{name}|{key}'] = wave_info(wp)
            tp = d / 'Calibrations' / f'Tilts_{key}.fits'
            if tp.exists():
                c['tilts'] = tilt_info(tp)
            fp = d / 'Calibrations' / f'Flat_{key}.fits'
            if fp.exists():
                with warnings.catch_warnings():
                    warnings.simplefilter('ignore')
                    c['flat'] = flat_info(fp, sp)
            r[key] = c
        res[name] = dict(calibs=r, science=science_info(d), boxes=sorted(box_slits(name)))
    with open(RESULTS / 'calibs.json', 'w') as fh:
        json.dump(res, fh, indent=1)
    hist = {}
    for key, hs in HISTS.items():
        if hs:
            h = np.sum(hs, axis=0).astype(float)
            h /= h.sum() * np.diff(HIST_EDGES)
            hist[key] = h
            hist[f'{key}_sigma'] = np.median([s['noise'][f'sigma_{key}'] for ds in res.values()
                                              for s in ds['science'] if 'science' in s['frametype']])
    if hist:
        np.savez(RESULTS / 'noise_hist.npz', edges=HIST_EDGES, **hist)
    with open(RESULTS / 'wavecal_lines.json', 'w') as fh:
        json.dump(all_lines, fh)
    print(json.dumps(res, indent=1)[:3000])


if __name__ == '__main__':
    main()
