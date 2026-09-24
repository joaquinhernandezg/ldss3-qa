#!/usr/bin/env python
"""
Build a QA bundle for a Magellan/LDSS3 reduction done with PypeIt.

Run this after ``run_pypeit`` finishes.  Point it at the reduction directory —
the one holding ``Calibrations/``, ``Science/`` and the ``.pypeit`` file — and it
writes a self-contained folder of plots plus a ``SUMMARY.txt`` of the numbers
that matter, and optionally a single ``.tar.gz`` to send back.

    ldss3_qa.py /path/to/magellan_ldss3_A
    ldss3_qa.py /path/to/magellan_ldss3_A -o ~/qa_2024run --tar

LDSS3-C reads one CCD through two amplifiers and writes each to its own file, so
several of these checks exist to confirm the two halves were joined correctly.
Nothing here modifies the reduction.
"""
import argparse
import re
import shutil
import sys
import tarfile
import traceback
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import AsinhNorm
from matplotlib.patches import Rectangle

SEAM = 1024          # amplifier boundary, trimmed unbinned frame
C1, C2, CG = '#2E6F9E', '#B0543C', '#B4881C'

notes = []
HOLES = set()        # spat_ids of apertures identified as alignment holes


def note(msg):
    notes.append(msg)
    print(f'    {msg}')


# --------------------------------------------------------------------- helpers
def find_setups(redux):
    """Find every (calib_key, label) in a reduction directory."""
    cal = redux / 'Calibrations'
    if not cal.is_dir():
        return []
    keys = sorted({m.group(1) for m in
                   (re.match(r'Slits_(.+)\.fits\.gz$', p.name) for p in cal.glob('Slits_*'))
                   if m})
    out = []
    for k in keys:
        grating = 'unknown'
        try:
            from pypeit.spec2dobj import AllSpec2DObj
            f2 = sorted((redux / 'Science').glob('spec2d_*.fits'))
            if f2:
                from astropy.io import fits
                hdr = fits.getheader(str(f2[0]), 0)
                grating = hdr.get('DISPNAME', hdr.get('GRISM', 'unknown'))
        except Exception:
            pass
        out.append((k, grating))
    return out


def load_spec2d(redux):
    from pypeit.spec2dobj import AllSpec2DObj
    from pypeit.specobjs import SpecObjs
    f2 = sorted((redux / 'Science').glob('spec2d_*.fits'))
    if not f2:
        return None, None, None
    # Prefer the longest exposure with a spec1d file: standards are reduced
    # without a global sky model and are not useful for the sky/noise checks.
    from astropy.io import fits
    have1d = [f for f in f2 if f.with_name(f.name.replace('spec2d_', 'spec1d_')).exists()] or f2
    pick = max(have1d, key=lambda f: float(fits.getheader(str(f), 0).get('EXPTIME', 0)))
    a = AllSpec2DObj.from_fits(str(pick), chk_version=False)
    s2 = a[a.detectors[0]]
    f1 = pick.with_name(pick.name.replace('spec2d_', 'spec1d_'))
    so = SpecObjs.from_fitsfile(str(f1), chk_version=False) if f1.exists() else None
    return s2, so, pick.name


def edges_of(slits):
    left, right, _ = slits.select_edges()
    return left, right


def save(fig, outdir, name):
    outdir.mkdir(parents=True, exist_ok=True)
    p = outdir / f'{name}.png'
    fig.savefig(p, dpi=110, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f'    wrote {name}.png')


def step_across(prof, seam=SEAM, gap=8, win=110):
    """Fit the trend either side of the seam and extrapolate, so a slope is not
    mistaken for a discontinuity.  Returns (percent step, left, right)."""
    xl = np.arange(max(seam - gap - win, 0), seam - gap)
    xr = np.arange(seam + gap, min(seam + gap + win, prof.size))
    ml = np.isfinite(prof[xl]) & (prof[xl] != 0)
    mr = np.isfinite(prof[xr]) & (prof[xr] != 0)
    if ml.sum() < 20 or mr.sum() < 20:
        return None, None, None
    a = np.polyfit(xl[ml], prof[xl][ml], 1)
    b = np.polyfit(xr[mr], prof[xr][mr], 1)
    el, er = np.polyval(a, seam), np.polyval(b, seam)
    return 100.0 * (er - el) / (0.5 * (el + er)), el, er


# ----------------------------------------------------------------------- plots
def plot_seam(redux, key, outdir):
    """Does the flux match where the two amplifier halves meet?"""
    from pypeit.flatfield import FlatImages
    from pypeit.slittrace import SlitTraceSet
    fp = redux / 'Calibrations' / f'Flat_{key}.fits'
    if not fp.exists():
        note('no Flat file, skipped the seam check')
        return
    fi = FlatImages.from_file(str(fp), chk_version=False)
    raw = fi.pixelflat_raw
    if raw is None:
        note('flat has no raw image, skipped the seam check')
        return
    prof = np.median(raw, axis=0)

    sp = redux / 'Calibrations' / f'Slits_{key}.fits.gz'
    spans = False
    if sp.exists():
        st = SlitTraceSet.from_file(str(sp), chk_version=False)
        mid = st.left_init.shape[0] // 2
        room = np.minimum(SEAM - st.left_init[mid], st.right_init[mid] - SEAM)
        spans = bool((room > 0).any())

    fig = plt.figure(figsize=(12.2, 4.4))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.9, 1.15], hspace=0.62, wspace=0.22)
    axf = fig.add_subplot(gs[:, 0])
    axr = fig.add_subplot(gs[0, 1])
    axn = fig.add_subplot(gs[1, 1])

    axf.plot(prof, lw=0.7, color='0.3')
    axf.axvline(SEAM, color=CG, lw=1.5, ls='--')
    axf.axvspan(0, SEAM, color=C1, alpha=0.07)
    axf.axvspan(SEAM, prof.size, color=C2, alpha=0.07)
    axf.set_xlim(0, prof.size)
    axf.set_xlabel('spatial pixel (trimmed)')
    axf.set_ylabel('flat, column median')
    axf.set_title('Spatial profile of the combined flat', fontsize=10)

    # The step is measured inside the slit that spans the seam: the fit windows
    # must not reach the slit edges, or the edge of the illumination is taken
    # for a discontinuity.
    gap, win = 8, 110
    if spans:
        win = int(min(win, np.max(room) - gap - 3))
        if win < 12:
            spans = False
            win = 110
    sl = slice(SEAM - gap - win, SEAM + gap + win)
    axr.plot(np.arange(SEAM - gap - win, SEAM + gap + win), prof[sl], lw=1.0, color='0.25')
    axr.axvline(SEAM, color=CG, lw=1.5, ls='--')

    jump = step_across(prof, win=win)[0] if spans else None
    if jump is None:
        axr.set_title('No slit spans enough of the seam', fontsize=10)
        axn.axis('off')
        note('no slit spans enough of the amplifier seam, so the flux step is not '
             'measurable')
    else:
        jump, el, er = step_across(prof, win=win)
        axr.plot([SEAM, SEAM], [el, er], color=CG, lw=2.6, solid_capstyle='butt')
        axr.set_title(f'Combined flat: {jump:+.2f} %', fontsize=10)
        pfn = getattr(fi, 'pixelflat_norm', None)
        jn = None
        if pfn is not None:
            jn, eln, ern = step_across(np.median(pfn, axis=0), win=win)
            if jn is not None:
                axn.plot(np.arange(SEAM - gap - win, SEAM + gap + win),
                         np.median(pfn, axis=0)[sl], lw=1.0, color='0.25')
                axn.axvline(SEAM, color=CG, lw=1.4, ls='--')
                axn.plot([SEAM, SEAM], [eln, ern], color=CG, lw=2.6, solid_capstyle='butt')
                axn.set_title(f'After normalising: {jn:+.2f} %', fontsize=10)
                axn.set_xlabel('spatial pixel (trimmed)')
        if jn is None:
            axn.axis('off')
        txt = f'seam step: raw flat {jump:+.2f} %'
        if jn is not None:
            txt += f'   normalised {jn:+.2f} %'
        axf.text(0.02, 0.03, txt, transform=axf.transAxes, fontsize=8.6,
                 family='monospace',
                 bbox=dict(fc='white', ec='0.8', boxstyle='round,pad=0.4'))
        note(f'amplifier seam: flux step {jump:+.2f} % in the combined flat'
             + (f', {jn:+.2f} % after normalising' if jn is not None else ''))
    fig.suptitle('Amplifier join continuity', fontsize=11.5, y=1.02)
    save(fig, outdir, '01_seam_continuity')


def plot_slits(redux, key, outdir):
    from pypeit.slittrace import SlitTraceSet
    sp = redux / 'Calibrations' / f'Slits_{key}.fits.gz'
    if not sp.exists():
        note('no Slits file, skipped the slit map')
        return
    st = SlitTraceSet.from_file(str(sp), chk_version=False)
    mid = st.left_init.shape[0] // 2
    nspec = st.left_init.shape[0]
    y = np.arange(nspec)
    fig, ax = plt.subplots(figsize=(11.5, 4.3))
    crossing = []
    for i in range(st.nslits):
        l, r = st.left_init[:, i], st.right_init[:, i]
        sp_ = l[mid] < SEAM < r[mid]
        if sp_:
            crossing.append(int(st.spat_id[i]))
        col = CG if sp_ else '0.35'
        ax.fill_betweenx(y, l, r, color=col, alpha=0.65 if sp_ else 0.3, lw=0)
        ax.text(0.5 * (l[mid] + r[mid]), nspec * 1.01, str(int(st.spat_id[i])),
                ha='center', fontsize=7, color=col, rotation=90)
    ax.axvline(SEAM, color=C2, lw=1.6, ls='--')
    ax.text(SEAM, nspec * 0.02, ' amplifier seam', color=C2, fontsize=9,
            fontweight='bold',
            bbox=dict(fc='white', ec='none', alpha=0.8, boxstyle='round,pad=0.2'))
    ax.set_xlim(0, st.nspat)
    ax.set_ylim(0, nspec)
    ax.set_xlabel('spatial pixel (trimmed)')
    ax.set_ylabel('spectral pixel')
    ax.set_title(f'{st.nslits} slits; {len(crossing)} cross the seam (gold)', fontsize=10.5)
    save(fig, outdir, '02_slit_map')
    nbad = int((st.mask != 0).sum())
    note(f'{st.nslits} slits traced, {nbad} flagged bad; seam-crossing slit ids {crossing}')


def plot_wavecal(redux, key, outdir):
    """Dispersion, line count and RMS per slit.  Dispersion is the discriminator:
    a slit whose dispersion differs from the rest has a wrong solution."""
    from pypeit.wavecalib import WaveCalib
    wp = redux / 'Calibrations' / f'WaveCalib_{key}.fits'
    if not wp.exists():
        note('no WaveCalib file, skipped the wavelength check')
        return
    wc = WaveCalib.from_file(str(wp), chk_version=False)
    sid, dwave, nlin, rms, wmin, wmax, fwhm = [], [], [], [], [], [], []
    for i, f in enumerate(wc.wv_fits):
        if f is None or f.wave_fit is None or f.wave_soln is None:
            continue
        ws = f.wave_soln
        sid.append(int(wc.spat_ids[i]))
        dwave.append(abs(ws[-1] - ws[0]) / ws.size)
        nlin.append(int(f.wave_fit.size))
        rms.append(float(f.rms))
        wmin.append(float(ws.min()))
        wmax.append(float(ws.max()))
        fwhm.append(float(f.fwhm))
    if not sid:
        note('no wavelength solutions found')
        return
    sid, dwave, nlin = np.array(sid), np.array(dwave), np.array(nlin)
    rms, fwhm = np.array(rms), np.array(fwhm)
    med = np.median(dwave)
    good = np.abs(dwave - med) / med < 0.05
    # Slitmask alignment holes are traced as slits.  They are much wider than a
    # slit in the dispersion direction, so their arc lines are broad.  Compare
    # with the narrow slits (a mask can have as many holes as slits, so the
    # median is not a safe reference).
    hole = fwhm > 1.5 * np.percentile(fwhm, 20)
    HOLES.update(int(x) for x in sid[hole])

    fig, axes = plt.subplots(1, 3, figsize=(12.2, 3.6))
    x = np.arange(sid.size)
    for ax, v, lab, t in ((axes[0], dwave, 'dispersion (Å/px)', 'Dispersion per slit'),
                          (axes[1], nlin, 'arc lines identified', 'Lines identified'),
                          (axes[2], rms, 'fit RMS (px)', 'Fit RMS')):
        ax.bar(x[good], v[good], color=C1, label='consistent')
        if (~good).any():
            ax.bar(x[~good], v[~good], color=C2, label='OUTLIER')
        ax.set_xticks(x)
        ax.set_xticklabels([str(s) for s in sid], rotation=90, fontsize=7)
        ax.set_xlabel('slit spat_id')
        ax.set_ylabel(lab)
        ax.set_title(t, fontsize=10)
        if ax is axes[0]:
            ax.axhline(med, color=CG, lw=1.2, ls='--')
            ax.legend(fontsize=7.5, frameon=False)
    fig.suptitle('Wavelength solution quality, slit by slit', fontsize=11.5, y=1.04)
    save(fig, outdir, '03_wavecal_quality')

    with open(outdir / '03_wavecal_table.txt', 'w') as fh:
        fh.write(f'{"spat_id":>8s} {"dWave":>8s} {"Nlin":>5s} {"fwhm":>6s} '
                 f'{"RMS(px)":>8s} {"wave_min":>9s} {"wave_max":>9s}  flag\n')
        for i in range(sid.size):
            flag = 'HOLE?' if hole[i] else ('ok' if good[i] else 'OUTLIER')
            fh.write(f'{sid[i]:8d} {dwave[i]:8.3f} {nlin[i]:5d} {fwhm[i]:6.2f} '
                     f'{rms[i]:8.3f} {wmin[i]:9.1f} {wmax[i]:9.1f}  {flag}\n')
    note(f'{int(good.sum())}/{sid.size} slits have a dispersion within 5% of the '
         f'median ({med:.3f} A/px); median fit RMS {np.median(rms):.3f} px')
    if hole.any():
        note(f'  apertures with broad arc lines, most likely alignment holes (ignore '
             f'them): {list(map(int, sid[hole]))}')
    if (~good & ~hole).any():
        note(f'  slits with a SUSPECT wavelength solution: {list(map(int, sid[~good & ~hole]))}')


def plot_spec2d(redux, outdir):
    s2, so, fname = load_spec2d(redux)
    if s2 is None:
        note('no spec2d files, skipped the 2D checks')
        return
    img = s2.sciimg - s2.skymodel
    gd = np.isfinite(img) & (s2.waveimg > 0)
    vals = img[gd] if gd.sum() > 1000 else img[np.isfinite(img)]
    lo, hi = np.percentile(vals, [3, 99.3])
    left, right = edges_of(s2.slits)
    y = np.arange(img.shape[0])
    mid = img.shape[0] // 2

    fig, ax = plt.subplots(figsize=(13.5, 6.4))
    ax.imshow(img, origin='lower', aspect='auto', cmap='gray_r',
              norm=AsinhNorm(linear_width=max((hi - lo) * 0.05, 1e-3), vmin=lo, vmax=hi),
              extent=[0, img.shape[1], 0, img.shape[0]])
    for i in range(s2.slits.nslits):
        sp_ = left[mid, i] < SEAM < right[mid, i]
        col = CG if sp_ else C1
        ax.plot(left[:, i], y, color=col, lw=0.9)
        ax.plot(right[:, i], y, color=col, lw=0.9)
        ax.text(0.5 * (left[mid, i] + right[mid, i]), img.shape[0] * 0.995,
                str(int(s2.slits.spat_id[i])), ha='center', va='top', fontsize=6.8,
                color=col, rotation=90,
                bbox=dict(fc='white', ec='none', alpha=0.72, boxstyle='round,pad=0.15'))
    nobj = 0
    if so is not None:
        for o in so:
            t = getattr(o, 'TRACE_SPAT', None)
            if t is None:
                continue
            nobj += 1
            ax.plot(t, np.arange(np.size(t)), color=C2, lw=1.1, ls='--')
    ax.axvline(SEAM, color=CG, lw=1.6, ls=':')
    ax.set_xlim(0, img.shape[1])
    ax.set_ylim(0, img.shape[0])
    ax.set_xlabel('spatial pixel  (amp 1 | amp 2)')
    ax.set_ylabel('spectral pixel')
    ax.set_title(f'Sky-subtracted 2D — {fname}\n{s2.slits.nslits} slit traces '
                 f'(gold crosses the seam), {nobj} object traces (dashed red)',
                 fontsize=10.5, pad=10)
    save(fig, outdir, '04_spec2d_overview')
    note(f'2D frame checked: {fname}, {nobj} objects extracted')

    # wavelength continuity across the seam
    w = s2.waveimg.copy()
    w[w <= 0] = np.nan
    steps = []
    for r in (np.linspace(0.12, 0.88, 5) * w.shape[0]).astype(int):
        j, _, _ = step_across(w[r], gap=6, win=90)
        if j is not None:
            row = w[r]
            _, el, er = step_across(row, gap=6, win=90)
            steps.append(er - el)
    fig, ax = plt.subplots(figsize=(7.6, 4.2))
    for r in (np.linspace(0.12, 0.88, 5) * w.shape[0]).astype(int):
        row = w[r]
        m = np.isfinite(row)
        if m.sum() > 50:
            ax.plot(np.arange(row.size)[m], row[m], lw=0.9, label=f'row {r}')
    ax.axvline(SEAM, color=C2, lw=1.6, ls=':')
    ax.set_xlabel('spatial pixel')
    ax.set_ylabel('wavelength (Å)')
    ax.legend(fontsize=7.5, frameon=False, ncol=2)
    if steps:
        ax.set_title(f'Wavelength across the seam: median jump {np.median(steps):+.3f} Å',
                     fontsize=10)
        note(f'wavelength jump across the seam = {np.median(steps):+.4f} A')
    else:
        ax.set_title('Wavelength across the seam: no slit spans it', fontsize=10)
    save(fig, outdir, '05_wavelength_across_seam')


def plot_spec1d(redux, outdir):
    from pypeit.specobjs import SpecObjs
    f1 = sorted((redux / 'Science').glob('spec1d_*.fits'))
    if not f1:
        note('no spec1d files, nothing was extracted')
        return
    so = SpecObjs.from_fitsfile(str(f1[0]), chk_version=False)
    rows = []
    for o in so:
        w = o.OPT_WAVE if o.OPT_WAVE is not None else o.BOX_WAVE
        c = o.OPT_COUNTS if o.OPT_COUNTS is not None else o.BOX_COUNTS
        if w is None or c is None:
            continue
        g = w > 1
        if g.sum() < 100:
            continue
        rows.append((o.NAME, w[g], c[g], float(np.median(c[g]))))
    if not rows:
        note('spec1d present but no usable extractions')
        return
    rows.sort(key=lambda r: -r[3])
    keep = rows[:6]
    fig, axes = plt.subplots(len(keep), 1, figsize=(11.5, 1.55 * len(keep)), sharex=True)
    axes = np.atleast_1d(axes)
    for ax, (name, w, c, m) in zip(axes, keep):
        ax.plot(w, c, lw=0.45, color=C1)
        ax.set_ylabel('counts', fontsize=8)
        ax.text(0.005, 0.86, f'{name}   median {m:,.0f}', transform=ax.transAxes,
                fontsize=8, family='monospace')
        ax.set_ylim(np.percentile(c, 0.5), np.percentile(c, 99.7) * 1.15)
    axes[-1].set_xlabel('wavelength (Å, vacuum)')
    fig.suptitle(f'Extracted spectra — {f1[0].name}', fontsize=11, y=0.997)
    fig.tight_layout()
    save(fig, outdir, '06_spec1d')
    note(f'{len(rows)} objects extracted in {f1[0].name}; '
         f'brightest median {keep[0][3]:,.0f} counts')


def plot_saturation(redux, key, outdir):
    """Alignment holes on a slitmask saturate and give unusable wavelength
    solutions.  Flag any slit whose flat is close to the ADC ceiling."""
    from pypeit.flatfield import FlatImages
    from pypeit.slittrace import SlitTraceSet
    fp = redux / 'Calibrations' / f'Flat_{key}.fits'
    sp = redux / 'Calibrations' / f'Slits_{key}.fits.gz'
    if not (fp.exists() and sp.exists()):
        return
    fi = FlatImages.from_file(str(fp), chk_version=False)
    raw = fi.pixelflat_raw
    if raw is None:
        return
    st = SlitTraceSet.from_file(str(sp), chk_version=False)

    # The ADC is 16-bit, so the ceiling in electrons is 65535 x gain.  Read the
    # gain from the reduction rather than assuming it: it differs between the
    # two amplifiers and with the readout mode (Slow / Fast / Turbo).
    gain = None
    try:
        s2, _, _ = load_spec2d(redux)
        if s2 is not None:
            gain = np.atleast_1d(s2.detector['gain']).astype(float)
    except Exception:
        pass
    if gain is None or gain.size == 0:
        gain = np.array([1.65])
        note('could not read the gain from the reduction; assuming 1.65 e-/ADU '
             'for the saturation check')
    ceiling = 65535.0 * float(np.max(gain))
    mid = raw.shape[0] // 2
    sid, peak = [], []
    for i in range(st.nslits):
        l, r = max(int(st.left_init[mid, i]), 0), min(int(st.right_init[mid, i]), raw.shape[1])
        if r - l < 4:
            continue
        sub = raw[mid - 400:mid + 400, l + 2:r - 2]
        if sub.size == 0:
            continue
        sid.append(int(st.spat_id[i]))
        peak.append(float(np.nanpercentile(sub, 90)))
    if not sid:
        return
    sid, peak = np.array(sid), np.array(peak)
    # Within 15 per cent of the ADC ceiling.  A lower threshold also catches
    # well-exposed science slits.
    hot = peak > 0.85 * ceiling
    fig, ax = plt.subplots(figsize=(9.6, 3.6))
    ax.bar(np.arange(sid.size)[~hot], peak[~hot], color=C1, label='normal')
    if hot.any():
        ax.bar(np.arange(sid.size)[hot], peak[hot], color=C2, label='near saturation')
    ax.axhline(ceiling, color=CG, lw=1.3, ls='--')
    ax.text(0, ceiling, f' ADC ceiling (65535 ADU x gain {np.max(gain):.2f})',
            color=CG, fontsize=8, va='bottom')
    ax.set_xticks(np.arange(sid.size))
    ax.set_xticklabels([str(s) for s in sid], rotation=90, fontsize=7)
    ax.set_xlabel('slit spat_id')
    ax.set_ylabel('flat 90th percentile (e-)')
    ax.set_title('Flat level per slit — saturated slits are usually alignment holes',
                 fontsize=10)
    ax.legend(fontsize=8, frameon=False)
    save(fig, outdir, '07_slit_saturation')
    if hot.any():
        note(f'slits near saturation (likely alignment holes, do not use for '
             f'science): {list(map(int, sid[hot]))}')


def plot_noise_model(redux, outdir):
    """Is the variance model right?  chi = (data - sky - object) x sqrt(ivar)
    should scatter with sigma 1.  Only the columns where the sky was actually
    modelled test this; elsewhere the residual still holds the full sky."""
    s2, _, fname = load_spec2d(redux)
    if s2 is None:
        return
    chi = (s2.sciimg - s2.skymodel - s2.objmodel) * np.sqrt(np.clip(s2.ivarmodel, 0, None))
    try:
        gd = (s2.bpmmask.flagged() == False) & np.isfinite(chi) & (s2.waveimg > 0)  # noqa: E712
    except Exception:
        gd = np.isfinite(chi) & (s2.waveimg > 0)
    sky = np.zeros(chi.shape[1])
    for c in range(chi.shape[1]):
        if gd[:, c].sum() > 200:
            sky[c] = np.median(s2.skymodel[:, c][gd[:, c]])
    sub = np.where(sky > 5)[0]
    if sub.size < 5:
        note('no sky-subtracted columns found, skipped the noise-model check')
        return
    # Exclude object-dominated pixels (standards, alignment stars): their residuals
    # measure the object-profile model, not the noise.
    var = np.where(s2.ivarmodel > 0, 1.0 / np.where(s2.ivarmodel > 0, s2.ivarmodel, 1), np.inf)
    gd = gd & (np.abs(s2.objmodel) < 0.1 * np.sqrt(var))
    if HOLES:
        gd &= ~np.isin(s2.slits.slit_img(initial=True), list(HOLES))
    v = chi[:, sub][gd[:, sub]]
    sg = 1.4826 * np.median(np.abs(v - np.median(v)))

    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    ax.hist(v, bins=180, range=(-6, 6), histtype='step', color=C1, lw=1.5,
            density=True, label=f'measured, sigma = {sg:.3f}')
    xs = np.linspace(-6, 6, 300)
    ax.plot(xs, np.exp(-xs ** 2 / 2) / np.sqrt(2 * np.pi), color='0.35', lw=1.1,
            ls='--', label='unit normal')
    ax.set_xlabel(r'$\chi$ = (data $-$ sky $-$ object) $\times\ \sqrt{\rm ivar}$')
    ax.set_ylabel('density')
    ax.set_title('Noise model, over sky-subtracted columns', fontsize=10.5)
    ax.legend(fontsize=8.5, frameon=False)
    rn = np.atleast_1d(s2.detector['ronoise'])
    ax.text(0.02, 0.02, f'read noise in use: {", ".join(f"{r:.2f}" for r in rn)} e-',
            transform=ax.transAxes, fontsize=8.2, family='monospace')
    save(fig, outdir, '08_noise_model')
    verdict = ('good' if 0.9 <= sg <= 1.15 else
               'noise UNDERSTATED' if sg > 1.15 else 'noise OVERSTATED')
    note(f'noise model: sigma(chi) = {sg:.3f} over sky-subtracted columns '
         f'({verdict}); read noise in use '
         f'{", ".join(f"{r:.2f}" for r in rn)} e-')


# ------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('redux', type=Path,
                    help='PypeIt reduction directory (contains Calibrations/ and Science/)')
    ap.add_argument('-o', '--outdir', type=Path, default=None,
                    help='where to write the bundle (default: <redux>/LDSS3_QA)')
    ap.add_argument('--tar', action='store_true',
                    help='also write a .tar.gz of the bundle, ready to send')
    ap.add_argument('--no-pypeit-qa', action='store_true',
                    help="do not copy PypeIt's own QA folder (makes the bundle smaller)")
    args = ap.parse_args()

    redux = args.redux.expanduser().resolve()
    if not (redux / 'Calibrations').is_dir():
        sys.exit(f'error: {redux} does not look like a PypeIt reduction directory '
                 '(no Calibrations/ inside)')
    out = (args.outdir or redux / 'LDSS3_QA').expanduser().resolve()
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    setups = find_setups(redux)
    if not setups:
        sys.exit(f'error: no Slits_*.fits.gz found in {redux}/Calibrations')

    print(f'reduction : {redux}')
    print(f'output    : {out}')
    for key, grating in setups:
        print(f'\n== calibration set {key}  (grating: {grating})')
        note(f'--- {key}  grating={grating} ---')
        sub = out / key
        for fn in (plot_seam, plot_slits, plot_wavecal, plot_saturation):
            try:
                fn(redux, key, sub)
            except Exception as e:
                note(f'{fn.__name__} FAILED: {type(e).__name__}: {e}')
                traceback.print_exc()
    print('\n== 2D and 1D spectra')
    for fn in (plot_spec2d, plot_spec1d, plot_noise_model):
        try:
            fn(redux, out / 'science')
        except Exception as e:
            note(f'{fn.__name__} FAILED: {type(e).__name__}: {e}')
            traceback.print_exc()

    # PypeIt's own QA, and the inputs a reviewer needs
    if not args.no_pypeit_qa and (redux / 'QA').is_dir():
        shutil.copytree(redux / 'QA', out / 'pypeit_qa')
        n = sum(1 for p in (out / 'pypeit_qa').rglob('*') if p.is_file())
        print(f'\n    copied {n} files of PypeIt QA')
    for pat in ('*.pypeit', '*.calib', '*.par'):
        for p in redux.glob(pat):
            shutil.copy(p, out / p.name)
    for p in redux.glob('*.log'):
        shutil.copy(p, out / p.name)

    with open(out / 'SUMMARY.txt', 'w') as fh:
        fh.write('LDSS3 reduction QA\n')
        fh.write(f'generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC}\n')
        fh.write(f'reduction {redux}\n')
        fh.write('=' * 68 + '\n\n')
        for n in notes:
            fh.write(n + '\n')
        fh.write('\nWhat to look at first:\n')
        fh.write('  * any line above saying SUSPECT, OUTLIER or FAILED\n')
        fh.write('  * the seam step: under ~0.5 % after normalising is fine\n')
        fh.write('  * the wavelength jump across the seam: should be < 0.1 A\n')
        fh.write('  * slits flagged near saturation or with broad arc lines are alignment holes, not science\n')
        fh.write('  * sigma(chi) should be close to 1.  Much above it means the\n')
        fh.write('    variance model understates the noise -- check the read noise.\n')

    print('\n--- summary ---')
    print((out / 'SUMMARY.txt').read_text())

    if args.tar:
        tgz = out.with_suffix('.tar.gz')
        with tarfile.open(tgz, 'w:gz') as tf:
            tf.add(out, arcname=out.name)
        print(f'bundle: {tgz}  ({tgz.stat().st_size / 1e6:.1f} MB)')
        print('Send that file back.')


if __name__ == '__main__':
    main()
