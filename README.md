# ldss3-qa

Quality-assessment scripts for Magellan/LDSS3 reductions done with the
[LDSS3 fork of PypeIt](https://github.com/joaquinhernandezg/PypeIt/tree/ldss3_multigrating).

Run this after `run_pypeit` finishes. It makes a folder of plots and a plain-text
summary, and can pack the lot into one `.tar.gz` to send back for a second pair
of eyes. It never modifies the reduction.

## Why this exists

LDSS3-C reads a single CCD through two amplifiers and writes each to its own
file (`ccd0042c1.fits`, `ccd0042c2.fits`). The fork joins them inside the
pipeline. Several of these checks exist specifically to confirm that join came
out right, which standard PypeIt QA does not test.

## Install

Nothing to install beyond the PypeIt environment you already use for the
reduction — the scripts import PypeIt itself.

```bash
git clone https://github.com/joaquinhernandezg/ldss3-qa.git
cd ldss3-qa
```

Activate the same conda environment you reduced in, then run it.

## Use

```bash
conda activate ldss3            # the env you reduced in
python ldss3_qa.py /path/to/magellan_ldss3_A --tar
```

The argument is your reduction directory — the one containing `Calibrations/`
and `Science/`. Options:

| Flag | Effect |
|---|---|
| `-o DIR` | write elsewhere (default `<redux>/LDSS3_QA`) |
| `--tar` | also write `LDSS3_QA.tar.gz`, ready to send |
| `--no-pypeit-qa` | skip copying PypeIt's own QA folder (much smaller bundle) |

Send back the `.tar.gz`. With PypeIt's own QA included it is usually 50–150 MB;
with `--no-pypeit-qa` it is a few MB.

## What you get

```
LDSS3_QA/
  SUMMARY.txt                     read this first
  <calib_key>/
    01_seam_continuity.png        flux match where the two amplifiers join
    02_slit_map.png               traced slits, seam marked
    03_wavecal_quality.png        dispersion / line count / RMS per slit
    03_wavecal_table.txt          the same, as numbers
    07_slit_saturation.png        which slits are saturated alignment boxes
  science/
    04_spec2d_overview.png        sky-subtracted 2D with slit + object traces
    05_wavelength_across_seam.png wavelength continuity across the join
    06_spec1d.png                 extracted spectra
    08_noise_model.png            is the variance model right?
  pypeit_qa/                      PypeIt's own QA, copied verbatim
  *.pypeit *.calib *.par *.log    the inputs, so the run can be reproduced
```

## Reading the summary

`SUMMARY.txt` ends with a short checklist. In order of importance:

**Anything saying SUSPECT, OUTLIER or FAILED.** A slit is flagged SUSPECT when
its dispersion differs from the median of the other slits by more than 5 %.
Dispersion is the reliable discriminator — fit RMS is not, because a badly wrong
solution can still have a small RMS if it latched onto few lines.

**Slits flagged near saturation.** On a slitmask these are almost always
alignment boxes. They saturate in the flats and arcs, their arc lines merge, and
their wavelength solutions are unusable. Do not extract science from them. This
is expected, not a failure.

**The seam step.** The combined flat normally shows a 1–3 % flux step where the
two amplifier halves meet — that is a small residual in the relative gain of the
two amplifiers. After normalising it should fall below about 0.5 %. If it does
not, the amplifier join or the flat is wrong.

**The noise model.** `sigma(chi)` is the scatter of
(data − sky − object) × √ivar over the columns where the sky was actually
modelled. It should be close to 1. Above ~1.15 the variance model is
understating the noise and the error bars on your spectra are too small; the
read noise is the usual cause.

**The wavelength jump across the seam.** Should be well under 0.1 Å. Values seen
in testing were 0.003–0.018 Å. A jump of an Ångström or more means the two
halves were misaligned when they were joined.

If no slit happens to lie across the seam, the last two checks report
"not measurable" rather than a number. That is normal for some slitmasks.

## Caveats

- Only 1×1 binning has been tested. The seam is assumed at spatial pixel 1024,
  which is correct for unbinned data; for binned data that constant in
  `ldss3_qa.py` needs dividing by the spatial binning.
- The saturation check compares against the ADC ceiling, which is 65535 ADU ×
  gain in electrons. The gain is read from your reduction, so it follows your
  readout mode; if it cannot be read the script says so and assumes 1.65. It is
  a flag, not a measurement.
- The script reads the first `spec2d`/`spec1d` pair it finds for the 2D and 1D
  checks, not all of them. For a multi-frame run it is a spot check.

## Commissioning report

`report/` holds a technical report on the LDSS3-C support in PypeIt, written
for LCO staff and users: what each reduction stage does, detector measurements
(read noise, gain, bad columns, dark current), wavelength-calibration
performance for the three grisms, example reductions, and recommendations.
The built PDF is `report/ldss3_pypeit_report.pdf`. Every number, table and
figure is produced by the scripts in `report/scripts/`; see the report's
appendix D, or run `make` in `report/`.
