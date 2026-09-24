"""
Build the PypeIt input files used for the commissioning report.

``pypeit_setup`` is run once over all of the raw data (see ``run_reductions.sh``).
This script takes the setup-specific files it writes and turns them into the
five datasets analysed in the report: it drops frames that are not used by the
default LDSS3 reduction (biases and darks, which are only used for the detector
measurements), fixes the few frames whose OBJECT card defeats the automatic
typing, groups calibrations by night and, for the multi-slit data, keeps a
representative subset of the science exposures.

Usage::

    python prepare_pypeit_files.py <setup_dir> <out_dir>

where ``<setup_dir>`` holds the ``magellan_ldss3_?`` directories written by
``pypeit_setup -c A,B,C,D,F,G`` and ``<out_dir>`` receives one directory per
dataset.
"""
import sys
from pathlib import Path

import numpy as np
from astropy.table import vstack

from pypeit.inputfiles import PypeItFile


def load(setup_dir, letter):
    return PypeItFile.from_file(setup_dir / f'magellan_ldss3_{letter}'
                                / f'magellan_ldss3_{letter}.pypeit')


def keep(tbl, names):
    """Rows whose filename stem (e.g. 'ccd0406') is in ``names``."""
    stems = np.array([f.split('c1.fits')[0].replace('ccd', '') for f in tbl['filename']])
    return tbl[np.isin(stems, names)]


def set_rows(tbl, names, column, value):
    stems = np.array([f.split('c1.fits')[0].replace('ccd', '') for f in tbl['filename']])
    tbl[column][np.isin(stems, names)] = value


def write(out_dir, name, pfile, tbl, cfg_lines=None):
    """Write a PypeIt file for one dataset."""
    config = {'rdx': {'spectrograph': 'magellan_ldss3'}}
    if cfg_lines is not None:
        config = cfg_lines
    tbl = tbl.copy()
    tbl.sort(['calib', 'frametype', 'filename'])
    tbl['calib'] = [str(c) for c in tbl['calib']]
    setup = {k: v for k, v in pfile.setup.items() if not k.startswith('Setup')}
    new = PypeItFile(config=config, file_paths=[str(p) for p in pfile.file_paths],
                     data_table=tbl, setup={f'Setup {name}': setup})
    d = out_dir / name
    d.mkdir(parents=True, exist_ok=True)
    new.write(d / f'{name}.pypeit')
    print(f'{name}: {len(tbl)} frames -> {d}')


def main(setup_dir, out_dir):
    setup_dir, out_dir = Path(setup_dir), Path(out_dir)

    # --- VPH-All, longslit, 2018 April (two nights) --------------------------
    p = load(setup_dir, 'A')
    t = p.data
    t = keep(t, ['0404', '0405', '0406', '0407', '0408',
                 '0422', '0423', '0424', '0425', '0450', '0451', '0452'])
    # "Hiltner600 - NeHeAr" is an arc; "nehear" is not one of the OBJECT
    # keywords recognised by the frame typing, so it is typed as a standard.
    set_rows(t, ['0424'], 'frametype', 'arc,tilt')
    t['calib'] = [0 if m < 58211.5 else 1 for m in t['mjd']]
    write(out_dir, 'vph_all_longslit', p, t)

    # --- VPH-Blue, longslit, 2019 February (calibrations only) ---------------
    p = load(setup_dir, 'B')
    t = p.data
    t = keep(t, ['2001', '2014', '2015', '2016', '2345', '1115', '1116', '1117'])
    t['calib'] = [0 if f[3:7] in ('2001', '2014', '2015', '2016') else 1
                  for f in t['filename']]
    write(out_dir, 'vph_blue_longslit', p, t)

    # --- VPH-Red, longslit, 2019 February (GD71 + LTT3218) --------------------
    p = load(setup_dir, 'C')
    t = p.data
    t = keep(t, ['1934', '1935', '1936', '1937', '1938', '1939', '1940', '1941',
                 '2004', '2005', '2006', '2007'])
    t['calib'] = [1 if f[3:7] in ('2004', '2005', '2006', '2007') else 0
                  for f in t['filename']]
    write(out_dir, 'vph_red_longslit', p, t)

    # --- VPH-Blue, multi-slit mask ATM3a2_1, 2022 November -------------------
    p = load(setup_dir, 'D')
    t = p.data
    t = keep(t, ['0024', '0025', '0026', '0027', '0028', '0029', '0045',
                 '0046', '0047', '0044', '0050', '0051', '0052'])
    t['calib'] = 0
    write(out_dir, 'vph_blue_mos', p, t)

    # --- VPH-Red, multi-slit mask CDFS01, 2025 December -----------------------
    # The flats were taken without the OG-590 blocking filter used for the arcs
    # and science frames, so pypeit_setup (correctly) splits them into two
    # setups.  They are combined here by hand; see the report for the
    # consequences.
    pf = load(setup_dir, 'F')
    pg = load(setup_dir, 'G')
    tf = keep(pf.data, [f'{i:04d}' for i in range(2011, 2021)])
    tg = keep(pg.data, ['2106', '2109', '2110', '2113', '2114', '2115'])
    tf['filter1'] = 'OG-590'
    t = vstack([tf, tg])
    t['calib'] = 0
    write(out_dir, 'vph_red_mos', pg, t)


if __name__ == '__main__':
    main(*sys.argv[1:3])
