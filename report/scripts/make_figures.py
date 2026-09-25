"""
Make every figure in the commissioning report.

Inputs: the raw data (``style.RAW``), the reductions (``style.REDUX``), the
PypeIt package data (archived templates and line lists), and the measurement
files written by ``measure_detector.py`` and ``measure_calibs.py`` in
``results/``.

Usage::

    python make_figures.py [name ...]      # default: all figures
"""
import sys
import json
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from astropy.table import Table

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style                                                       # noqa: E402
from style import (FIGS, RESULTS, REDUX, RAW, DATASETS, NO_OBJECT, TEXTWIDTH, AMP_COLOR, ION_COLOR,
                   BLUE, ORANGE, AQUA, YELLOW, VIOLET, INK, INK2, GRID, save)  # noqa: E402

from pypeit import dataPaths                                       # noqa: E402
from pypeit.core.wavecal import waveio                             # noqa: E402

SEAM = 1024
GRISMS = ['VPH-Blue', 'VPH-Red', 'VPH-All']
TEMPLATES = {'VPH-Blue': 'magellan_ldss3_vph_blue_HeINeIArI.fits',
             'VPH-Red': 'magellan_ldss3_vph_red_HeINeIArI.fits',
             'VPH-All': 'magellan_ldss3_vph_all_HeINeIArI.fits'}
LONGSLIT = {'VPH-Blue': 'vph_blue_longslit', 'VPH-Red': 'vph_red_longslit',
            'VPH-All': 'vph_all_longslit'}


def load_json(name):
    with open(RESULTS / name) as fh:
        return json.load(fh)


def calib_dir(ds):
    return REDUX / ds / 'Calibrations'


def first(ds, prefix, suffix='.fits'):
    return sorted(calib_dir(ds).glob(f'{prefix}_*{suffix}'))[0]


def template(grism):
    path = dataPaths.reid_arxiv.get_file_path(TEMPLATES[grism])
    return Table.read(path)


# ---------------------------------------------------------------- detector ---
def fig_detector_layout():
    """Raw amplifier files and the assembled, trimmed frame."""
    from astropy.io import fits
    from pypeit.spectrographs.util import load_spectrograph
    spec = load_spectrograph('magellan_ldss3')
    f1 = RAW / 'e190219_20' / 'ccd1940c1.fits'
    f2 = RAW / 'e190219_20' / 'ccd1940c2.fits'
    d1, d2 = fits.getdata(f1).astype(float), fits.getdata(f2).astype(float)
    _, arr, _, _, dsec, osec = spec.get_rawimage(str(f1), 1)

    fig, axes = plt.subplots(1, 3, figsize=(TEXTWIDTH, 3.3),
                             gridspec_kw=dict(width_ratios=[1.1, 1.1, 2.4], wspace=0.55))
    norm = AsinhNorm(linear_width=300, vmin=np.percentile(d1, 1), vmax=np.percentile(d1, 99.7))
    for ax, d, lab in ((axes[0], d1, r'\texttt{ccd1940c1.fits}'), (axes[1], d2, r'\texttt{ccd1940c2.fits}')):
        ax.imshow(d, cmap='Greys', norm=norm, aspect='auto', interpolation='nearest',
                  rasterized=True)
        ax.set_title(lab.replace(r'\texttt{', '').replace('}', ''), fontsize=7.5)
        ax.set_xlabel('NAXIS1 (pix)')
    axes[0].set_ylabel('NAXIS2 (pix)')
    # Assembled frame: (spectral, spatial) with overscan regions outlined
    a = np.ma.masked_where(dsec == 0, arr)
    axes[2].imshow(arr.T, cmap='Greys', norm=norm, aspect='auto', interpolation='nearest',
                   rasterized=True)
    for amp in (1, 2):
        cols = np.where((osec == amp).any(axis=0))[0]
        axes[2].axhspan(cols.min(), cols.max(), color=AMP_COLOR[amp], alpha=0.25, lw=0)
        dcols = np.where((dsec == amp).any(axis=0))[0]
        axes[2].text(100, 0.5 * (dcols.min() + dcols.max()), f'amp {amp}', ha='left',
                     va='center', fontsize=7.5, color=AMP_COLOR[amp],
                     bbox=dict(fc='white', ec='none', alpha=0.8, pad=1))
    axes[2].axhline(1152, color=INK, lw=0.6, ls='--')
    axes[2].set_xlabel('spectral pixel')
    axes[2].set_ylabel('spatial pixel (with overscan)')
    axes[2].set_title('assembled by PypeIt (transposed)', fontsize=7.5)
    save(fig, 'detector_layout')


def fig_bias_bpm():
    det = load_json('detector.json')
    prof = np.load(RESULTS / 'bias_profiles.npz')
    keys = sorted(prof.files)
    fig, axes = plt.subplots(2, 1, figsize=(TEXTWIDTH, 3.6), sharey=False,
                             gridspec_kw=dict(height_ratios=[1, 1], hspace=0.35))
    x = np.arange(2048)
    pick = [k for k in keys if k.startswith('2018')][:1] + [k for k in keys if k.startswith('2022')][:1] \
        + [k for k in keys if k.startswith('2025')][:1]
    colors = [BLUE, ORANGE, AQUA]
    for ax, xlim in ((axes[0], (0, 2048)), (axes[1], (1380, 1720))):
        for (c1, c2) in det['bpm_ranges']:
            ax.axvspan(c1 - 0.5, c2 + 0.5, color=YELLOW, alpha=0.35, lw=0)
        for k, c in zip(pick, colors):
            p = prof[k]
            p = p - np.median(p[20:1004])
            ax.plot(x, p, color=c, lw=0.7, label=k.split(' ')[0])
        ax.axvline(SEAM - 0.5, color=INK2, lw=0.6, ls='--')
        ax.set_xlim(*xlim)
        ax.set_ylim(-6, 6)
        ax.set_ylabel('bias $-$ median (ADU)')
    axes[1].set_xlabel('spatial pixel (trimmed frame)')
    axes[0].legend(ncol=3, loc='lower center')
    axes[0].text(0.01, 0.9, 'shaded: masked columns', transform=axes[0].transAxes, fontsize=7,
                 color=INK2)
    save(fig, 'bias_bpm')


def fig_readnoise():
    det = load_json('detector.json')
    bias = det['bias']
    keys = sorted(bias)
    fig, ax = plt.subplots(figsize=(TEXTWIDTH * 0.62, 2.6))
    x = np.arange(len(keys))
    for amp, off in ((1, -0.08), (2, 0.08)):
        y = [bias[k][str(amp)]['ron_e'] for k in keys]
        e = [bias[k][str(amp)]['ron_e_std'] for k in keys]
        ax.errorbar(x + off, y, yerr=e, fmt='o', ms=4, color=AMP_COLOR[amp], capsize=0,
                    lw=0.8, label=f'amp {amp}, measured')
    enoise = sorted({(v[2], v[4]) for s in det['header_survey'].values() for v in s})[0]
    for amp, ls in ((1, '-'), (2, '--')):
        ax.axhline({1: 7.0, 2: 7.2}[amp], color=AMP_COLOR[amp], lw=0.7, ls=ls)
        ax.axhline(enoise[amp - 1], color=AMP_COLOR[amp], lw=0.7, ls=ls, alpha=0.45)
    ax.text(len(keys) - 0.45, 7.3, 'published, Fast mode', fontsize=6.5, ha='right', color=INK2)
    ax.text(len(keys) - 0.45, 5.15, 'ENOISE keyword', fontsize=6.5, ha='right', color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels([k.split(' ')[0] for k in keys], rotation=30, ha='right')
    ax.set_ylabel('read noise (e$^-$)')
    ax.set_xlim(-0.5, len(keys) - 0.5)
    ax.set_ylim(4.0, 8.0)
    ax.legend(loc='lower left', fontsize=6.5, ncol=2)
    save(fig, 'readnoise')


def fig_ptc():
    det = load_json('detector.json')
    ptc = np.load(RESULTS / 'ptc.npz')
    # Use the set with the most points on the fit
    # The accepted set with the most points (selected in make_tables.py)
    best = (RESULTS / 'ptc_best.txt').read_text().strip()
    fig, ax = plt.subplots(figsize=(TEXTWIDTH * 0.62, 2.6))
    for amp in (1, 2):
        S, V, good = ptc[f'{best}|amp{amp}']
        good = good.astype(bool)
        g = det['ptc'][best][str(amp)]['gain']
        ge = det['ptc'][best][str(amp)]['gain_err']
        ax.plot(S[good], V[good], 'o', ms=3, color=AMP_COLOR[amp],
                label=f'amp {amp}: $g = {g:.2f}\\pm{ge:.2f}$ e$^-$/ADU')
        ax.plot(S[~good], V[~good], 'o', ms=3, mfc='none', color=AMP_COLOR[amp], alpha=0.5)
        s = np.linspace(0, S[good].max() * 1.05, 10)
        ax.plot(s, s / g, color=AMP_COLOR[amp], lw=0.8)
    ax.set_xlabel('mean signal (ADU)')
    ax.set_ylabel(r'$\sigma^2_{\rm pair}/2 - \sigma^2_{\rm RN}$ (ADU$^2$)')
    ax.set_xlim(0, None)
    ax.set_ylim(0, None)
    ax.legend(loc='upper left', fontsize=6.5)
    save(fig, 'ptc')
    return best


# ------------------------------------------------------------ slit tracing ---
def fig_slit_edges():
    from pypeit.edgetrace import EdgeTraceSet
    from pypeit.slittrace import SlitTraceSet
    panels = [('vph_red_longslit', 'VPH-Red\nlongslit'),
              ('vph_blue_mos', 'VPH-Blue\nmask ATM3a2_1'),
              ('vph_red_mos', 'VPH-Red\nmask CDFS01')]
    cal = load_json('calibs.json')
    fig, axes = plt.subplots(1, 3, figsize=(TEXTWIDTH, 3.4), sharey=True,
                             gridspec_kw=dict(wspace=0.05))
    for ax, (ds, lab) in zip(axes, panels):
        boxes = set(cal[ds]['boxes'])
        et = EdgeTraceSet.from_file(str(first(ds, 'Edges', '.fits.gz')), chk_version=False)
        st = SlitTraceSet.from_file(str(first(ds, 'Slits', '.fits.gz')), chk_version=False)
        img = et.traceimg.image
        v1, v2 = np.percentile(img[img > 0], [2, 99.5])
        ax.imshow(img, cmap='Greys', norm=AsinhNorm(linear_width=v2 / 10, vmin=0, vmax=v2),
                  aspect='auto', interpolation='nearest', rasterized=True)
        y = np.arange(st.nspec)
        for i in range(st.nslits):
            c = ORANGE if int(st.spat_id[i]) in boxes else BLUE
            ax.plot(st.left_init[:, i], y, color=c, lw=0.6)
            ax.plot(st.right_init[:, i], y, color=c, lw=0.6, ls='--')
        ax.axvline(SEAM - 0.5, color=YELLOW, lw=0.7)
        nb = sum(int(x) in boxes for x in st.spat_id)
        ax.set_title(f'{lab}: {st.nslits - nb} slits' + (f' + {nb} holes' if nb else ''), fontsize=7)
        ax.set_xlabel('spatial pixel')
    axes[0].set_ylabel('spectral pixel')
    save(fig, 'slit_edges')


# ------------------------------------------------------------ flat field ---
def fig_flat():
    from pypeit.flatfield import FlatImages
    from pypeit.slittrace import SlitTraceSet
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
    from ldss3_qa import step_across
    ds = 'vph_red_longslit'
    fi = FlatImages.from_file(str(first(ds, 'Flat')), chk_version=False)
    fig = plt.figure(figsize=(TEXTWIDTH, 3.6))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.2, 1.2, 1.3], wspace=0.4, hspace=0.45)
    ax0 = fig.add_subplot(gs[:, 0])
    ax1 = fig.add_subplot(gs[:, 1], sharey=ax0)
    ax2 = fig.add_subplot(gs[0, 2])
    ax3 = fig.add_subplot(gs[1, 2])
    raw = fi.pixelflat_raw
    ax0.imshow(raw, cmap='Greys_r', aspect='auto', interpolation='nearest', rasterized=True,
               vmin=0, vmax=np.percentile(raw, 99.5))
    ax0.set_title('combined flat (e$^-$)', fontsize=7.5)
    norm = fi.pixelflat_norm
    ax1.imshow(norm, cmap='RdBu_r', vmin=0.97, vmax=1.03, aspect='auto',
               interpolation='nearest', rasterized=True)
    ax1.set_title('pixel flat, 0.97--1.03'.replace('--', '–'), fontsize=7.5)
    for ax in (ax0, ax1):
        ax.set_xlabel('spatial pixel')
    ax0.set_ylabel('spectral pixel')
    mid = raw.shape[0] // 2
    for ax, img, lab in ((ax2, raw, 'combined flat'), (ax3, norm, 'normalised pixel flat')):
        prof = np.median(img[mid - 1000:mid + 1000], axis=0)
        sl = slice(SEAM - 120, SEAM + 120)
        xs = np.arange(prof.size)[sl]
        p = prof[sl] / np.median(prof[sl])
        ax.plot(xs, p, color=INK, lw=0.8)
        ax.axvline(SEAM - 0.5, color=YELLOW, lw=0.8)
        cal = load_json('calibs.json')
        f = list(cal[ds]['calibs'].values())[0]['flat']
        st_ = SlitTraceSet.from_file(str(first(ds, 'Slits', '.fits.gz')), chk_version=False)
        room = max(min(SEAM - st_.left_init[mid, i], st_.right_init[mid, i] - SEAM)
                   for i in range(st_.nslits))
        j, _, _ = step_across(prof, win=min(int(room) - 11, 110))
        ax.set_title(f'{lab}: step {j:+.2f}\\%'.replace('\\%', '%'), fontsize=7.5)
        ax.set_ylabel('relative level')
    ax3.set_xlabel('spatial pixel')
    save(fig, 'flat')


# ---------------------------------------------------- wavelength templates ---
def _identified_lines(ds):
    lines = load_json('wavecal_lines.json')
    key = [k for k in lines if k.startswith(ds + '|')][0]
    rows = lines[key]
    # choose the slit with the most lines
    by = {}
    for r in rows:
        by.setdefault(r['spat_id'], []).append(r)
    sid = max(by, key=lambda s: sum(r['used'] for r in by[s]))
    return sid, by[sid]


def fig_arc_templates():
    from matplotlib.transforms import blended_transform_factory
    fig, axes = plt.subplots(3, 1, figsize=(TEXTWIDTH, 7.8), gridspec_kw=dict(hspace=0.55))
    for ax, grism in zip(axes, GRISMS):
        t = template(grism)
        w, f = np.asarray(t['wave']), np.asarray(t['flux'])
        f = f / np.nanmax(f)
        ax.plot(w, np.clip(f, 1e-4, None), color=INK, lw=0.5)
        ax.set_yscale('log')
        ax.set_ylim(3e-4, 1.2)
        tr = blended_transform_factory(ax.transData, ax.transAxes)
        _, rows = _identified_lines(LONGSLIT[grism])
        used = [r for r in rows if r['used']]
        for r in used:
            ion = r['ion'].replace(' ', '')
            ax.plot([r['wave']] * 2, [1.02, 1.09], color=ION_COLOR.get(ion, INK2), lw=0.8,
                    transform=tr, clip_on=False)
        # label the brightest lines, keeping labels apart
        amp = np.interp([r['wave'] for r in used], w, f)
        order = np.argsort(amp)[::-1]
        placed = []
        span = w.max() - w.min()
        for i in order:
            x = used[i]['wave']
            if any(abs(x - p) < 0.035 * span for p in placed):
                continue
            placed.append(x)
            ax.text(x, 1.11, f"{used[i]['ion'].replace(' ', '')} {x:.1f}", rotation=90,
                    fontsize=5.2, ha='center', va='bottom', color=INK2, transform=tr)
            if len(placed) >= 14:
                break
        ax.set_xlim(w.min(), w.max())
        ax.set_ylabel('normalised counts')
        ax.text(0.005, 0.95, f'{grism}: {len(used)} lines in the fit', transform=ax.transAxes,
                fontsize=7.5, va='top')
    axes[-1].set_xlabel('vacuum wavelength (Å)')
    handles = [plt.Line2D([], [], color=ION_COLOR[i], lw=1.2, label=i) for i in ('HeI', 'NeI', 'ArI')]
    axes[-1].legend(handles=handles, ncol=3, loc='lower right', fontsize=6.5)
    save(fig, 'arc_templates')


def fig_wavecal_residuals():
    fig, axes = plt.subplots(3, 1, figsize=(TEXTWIDTH, 4.6), sharey=True,
                             gridspec_kw=dict(hspace=0.5))
    for ax, grism in zip(axes, GRISMS):
        sid, rows = _identified_lines(LONGSLIT[grism])
        lim = 1.2
        for r in rows:
            ion = r['ion'].replace(' ', '')
            y = np.clip(r['resid'], -0.95 * lim, 0.95 * lim)
            mk = ('o' if r['used'] else 'x') if abs(r['resid']) < lim else ('^' if r['resid'] > 0 else 'v')
            ax.plot(r['wave'], y, mk, ms=3,
                    color=ION_COLOR.get(ion, INK2), mfc=ION_COLOR.get(ion, INK2) if r['used'] else 'none')
        ax.set_ylim(-lim, lim)
        res = np.array([r['resid'] for r in rows if r['used']])
        rms = np.sqrt(np.mean(res ** 2))
        ax.axhline(0, color=INK2, lw=0.5)
        ax.axhspan(-rms, rms, color=GRID, alpha=0.5, lw=0)
        ax.set_title(f'{grism}, slit {sid}: {res.size} lines, RMS {rms:.3f} Å', fontsize=7.5)
        ax.set_ylabel('residual (Å)')
    axes[-1].set_xlabel('vacuum wavelength (Å)')
    save(fig, 'wavecal_residuals')


def fig_sky_template():
    t = Table.read(dataPaths.reid_arxiv.get_file_path('magellan_ldss3_vph_red_sky.fits'))
    oh = waveio.load_line_lists(['OH_LDSS3_vac'])[0]
    w, f = np.asarray(t['wave']), np.asarray(t['flux'])
    f = f / np.nanpercentile(f, 99.5)
    fig, ax = plt.subplots(figsize=(TEXTWIDTH, 2.3))
    ax.plot(w, f, color=INK, lw=0.5)
    sel = (oh['wave'] > w.min()) & (oh['wave'] < w.max())
    from matplotlib.transforms import blended_transform_factory
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    for x in oh['wave'][sel]:
        ax.plot([x, x], [1.02, 1.1], color=VIOLET, lw=0.5, transform=tr, clip_on=False)
    ax.set_xlim(w.min(), w.max())
    ax.set_ylim(-0.02, 1.05)
    ax.set_xlabel('vacuum wavelength (Å)')
    ax.set_ylabel('normalised counts')
    ax.text(0.005, 0.95, f'VPH-Red sky template; ticks: {sel.sum()} lines of OH_LDSS3_vac',
            transform=ax.transAxes, fontsize=7.5, va='top')
    save(fig, 'sky_template')


def fig_mos_wavecal():
    cal = load_json('calibs.json')
    fig, axes = plt.subplots(2, 3, figsize=(TEXTWIDTH, 4.0), gridspec_kw=dict(hspace=0.6, wspace=0.45))
    for row, ds in zip(axes, ('vph_blue_mos', 'vph_red_mos')):
        c = list(cal[ds]['calibs'].values())[0]
        boxes = set(cal[ds]['boxes'])
        sl = c['wave']['slits']
        x = np.arange(len(sl))
        isbox = np.array([s['spat_id'] in boxes for s in sl])
        for ax, key, lab in zip(row, ('disp', 'nlin', 'rms_ang'),
                                ('dispersion (Å/pix)', 'lines in fit', 'RMS (Å)')):
            v = np.array([s.get(key, np.nan) if s.get('ok') else np.nan for s in sl])
            ax.bar(x[~isbox], v[~isbox], color=BLUE, width=0.8, label='slit')
            ax.bar(x[isbox], v[isbox], color=ORANGE, width=0.8, label='alignment hole')
            for xi in x[isbox & ~np.isfinite(v)]:
                ax.text(xi, 0.02, 'x', ha='center', va='bottom', color=ORANGE, fontsize=7,
                        transform=ax.get_xaxis_transform())
            ax.set_ylabel(lab)
            ax.set_xticks(x)
            ax.set_xticklabels([str(s['spat_id']) for s in sl], rotation=90, fontsize=5.5)
        med = np.median([s['disp'] for s, b in zip(sl, isbox) if s.get('ok') and not b])
        row[0].set_ylim(med * 0.9, med * 1.1)
        row[2].set_ylim(0, 0.4)
        row[0].set_title(DATASETS[ds][2], fontsize=7.5, loc='left')
    axes[1, 1].set_xlabel('slit (spat\\_id)'.replace('\\', ''))
    axes[0, 2].legend(fontsize=6, loc='upper right')
    save(fig, 'mos_wavecal')


# -------------------------------------------------------- science examples ---
def fig_science_examples():
    """For each data set with science or standard frames: the sky-subtracted 2D
    spectrum of the slit holding the best object (holes excluded), and its
    extracted spectrum.  Faint multi-slit targets are shown binned, with the
    sky spectrum (scaled) for reference."""
    from pypeit.spec2dobj import AllSpec2DObj
    from pypeit.specobjs import SpecObjs
    cal = load_json('calibs.json')
    sets = [ds for ds in ('vph_all_longslit', 'vph_red_longslit', 'vph_blue_mos')
            if ds not in NO_OBJECT and cal.get(ds, {}).get('science')]
    fig, axes = plt.subplots(len(sets), 2, figsize=(TEXTWIDTH, 1.95 * len(sets)),
                             gridspec_kw=dict(width_ratios=[1, 2.1], hspace=0.65, wspace=0.32))
    for (ax2, ax1), ds in zip(np.atleast_2d(axes), sets):
        sci = cal[ds]['science']
        best = None
        for s in sci:
            for o in s.get('objects', []):
                if o.get('in_box') or o['sn_median'] is None:
                    continue
                if best is None or o['sn_median'] > best[1]['sn_median']:
                    best = (s['file'], o)
        f2 = REDUX / ds / 'Science' / best[0]
        a = AllSpec2DObj.from_fits(str(f2), chk_version=False)
        s2 = a[a.detectors[0]]
        so = SpecObjs.from_fitsfile(str(f2).replace('spec2d_', 'spec1d_'), chk_version=False)
        obj = [o for o in so if o.NAME == best[1]['name']][0]
        # 2D: the slit holding the object, rectified along the slit centre
        slitmask = s2.slits.slit_img(initial=True)
        cols = np.where((slitmask == obj.SLITID).any(axis=0))[0]
        c = int(round(obj.SPAT_PIXPOS))
        half = min(35, (cols[-1] - cols[0]) // 2 + 3)
        rows = np.arange(s2.sciimg.shape[0])
        trace = np.round(obj.TRACE_SPAT).astype(int)
        resid = s2.sciimg - s2.skymodel
        rect = np.array([resid[r, t - half:t + half] for r, t in zip(rows, trace)])
        wave = s2.waveimg[rows, trace]
        good = wave > 0
        vmax = np.percentile(rect[good], 99.5)
        ax2.imshow(rect[good].T, aspect='auto', cmap='Greys', vmin=-0.3 * vmax, vmax=vmax,
                   extent=[wave[good].min(), wave[good].max(), -half * 0.189, half * 0.189],
                   interpolation='nearest', rasterized=True, origin='lower')
        ax2.set_ylabel('offset (arcsec)')
        ax2.set_title(DATASETS[ds][2], fontsize=7, loc='left')
        m = obj.OPT_MASK & (obj.OPT_WAVE > 0)
        wv, fl = obj.OPT_WAVE[m], obj.OPT_COUNTS[m]
        sig = 1 / np.sqrt(obj.OPT_COUNTS_IVAR[m])
        faint = best[1]['sn_median'] < 5
        if faint:
            nb = 8
            n = wv.size // nb * nb
            wv = wv[:n].reshape(-1, nb).mean(1)
            iv = obj.OPT_COUNTS_IVAR[m][:n].reshape(-1, nb)
            fl = (obj.OPT_COUNTS[m][:n].reshape(-1, nb) * iv).sum(1) / iv.sum(1)
            sig = 1 / np.sqrt(iv.sum(1))
            sky = obj.OPT_COUNTS_SKY[m][:n].reshape(-1, nb).mean(1)
            scale = np.percentile(fl, 99) / np.percentile(sky, 99)
            ax1.plot(wv, sky * scale, color=GRID, lw=0.6, label='sky (scaled)')
        ax1.plot(wv, fl, color=BLUE, lw=0.6, label='counts' + (f' ({nb}-pix bins)' if faint else ''))
        ax1.plot(wv, sig, color=ORANGE, lw=0.6, label=r'$\sigma$')
        lo = min(0, np.percentile(fl, 1) * 1.2)
        ax1.set_ylim(lo, np.percentile(fl, 99.5) * 1.25)
        ax1.set_xlim(wv.min(), wv.max())
        ax1.set_title(f'slit {obj.SLITID}, median S/N per pixel {best[1]["sn_median"]:.1f}',
                      fontsize=7, loc='left')
        ax1.set_ylabel('counts (e$^-$)')
        ax1.legend(loc='upper right', fontsize=5.5, ncol=3)
    np.atleast_2d(axes)[-1, 0].set_xlabel('wavelength (Å)')
    np.atleast_2d(axes)[-1, 1].set_xlabel('wavelength (Å)')
    save(fig, 'science_examples')


def fig_noise_model():
    cal = load_json('calibs.json')
    hist = np.load(RESULTS / 'noise_hist.npz')
    fig, ax = plt.subplots(figsize=(TEXTWIDTH * 0.62, 2.6))
    edges = hist['edges']
    x = 0.5 * (edges[1:] + edges[:-1])
    for key, color, lab in (('published', BLUE, 'read noise from readout mode'),
                            ('enoise', ORANGE, 'read noise from ENOISE card')):
        h = hist[key]
        s = float(hist[f'{key}_sigma'])
        ax.step(x, h, where='mid', color=color, lw=0.9, label=f'{lab}: $\\sigma_\\chi = {s:.3f}$')
    ax.plot(x, np.exp(-x ** 2 / 2) / np.sqrt(2 * np.pi), color=INK2, lw=0.7, ls='--',
            label='unit normal')
    ax.set_yscale('log')
    ax.set_ylim(1e-4, 1)
    ax.set_xlabel(r'$\chi$ = (data $-$ sky $-$ object) $\sqrt{\mathrm{ivar}}$')
    ax.set_ylabel('density')
    ax.legend(fontsize=6, loc='upper left', bbox_to_anchor=(0.0, 1.0))
    ax.set_ylim(1e-4, 3)
    save(fig, 'noise_model')


ALL = {f.__name__[4:]: f for f in (fig_detector_layout, fig_bias_bpm, fig_readnoise, fig_ptc,
                                   fig_slit_edges, fig_flat, fig_arc_templates,
                                   fig_wavecal_residuals, fig_sky_template, fig_mos_wavecal,
                                   fig_science_examples, fig_noise_model)}

if __name__ == '__main__':
    style.setup()
    names = sys.argv[1:] or list(ALL)
    for n in names:
        try:
            ALL[n]()
        except Exception as e:                                     # noqa: BLE001
            print(f'FAILED {n}: {type(e).__name__}: {e}')
            if len(sys.argv) > 1:
                raise
