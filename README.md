# ldss3-qa

Quality assessment for Magellan/LDSS3 spectroscopic reductions made with the
[LDSS3 fork of PypeIt](https://github.com/joaquinhernandezg/PypeIt/tree/ldss3_multigrating).

Run one command after `run_pypeit` finishes. You get a folder of diagnostic
plots and a plain-text summary that says, in words, whether each check passed —
plus an optional single archive to send back for a second opinion.

It reads the reduction and writes elsewhere. **It never modifies your data.**

---

## Contents

- [Why this exists](#why-this-exists)
- [Before you start](#before-you-start)
- [Install](#install)
- [Check your setup](#check-your-setup)
- [Run it](#run-it)
- [What you get](#what-you-get)
- [Reading the results](#reading-the-results)
- [Troubleshooting](#troubleshooting)
- [Sending results back](#sending-results-back)
- [Commissioning report](#commissioning-report)
- [Caveats](#caveats)

---

## Why this exists

LDSS3-C reads a single CCD through two amplifiers and writes **each amplifier to
its own file** (`ccd0042c1.fits` and `ccd0042c2.fits`). The fork joins them
inside the pipeline, so there is no external merging step.

That join is the thing most worth checking, and standard PypeIt QA does not test
it. These scripts do: whether the flux matches where the two halves meet,
whether the wavelength solution runs continuously across the boundary, and
whether the noise model holds on both sides. They also catch the two problems
that turn up most often in LDSS3 slitmask data — slits with a wrong wavelength
solution, and alignment holes being reduced as if they were science slits.

---

## Before you start

You need two things:

1. **A completed PypeIt reduction** — a directory containing `Calibrations/` and
   `Science/`, usually named something like `magellan_ldss3_A`.
2. **The conda environment you reduced in.** These scripts import PypeIt itself,
   so there is nothing extra to install; but they must run somewhere the LDSS3
   fork is importable.

Not reduced anything yet? Start with
[`LDSS3_REDUCTION_MANUAL.md`](LDSS3_REDUCTION_MANUAL.md) in this repository,
which covers installing the fork and running a reduction from scratch.

---

## Install

```bash
git clone https://github.com/joaquinhernandezg/ldss3-qa.git
cd ldss3-qa
```

That is the whole installation — no package to build, no dependencies to
resolve. Everything needed comes from PypeIt.

To update later:

```bash
cd ldss3-qa && git pull
```

---

## Check your setup

Before committing to a long reduction, confirm the environment is right:

```bash
conda activate ldss3          # the environment you reduce in
python ldss3_qa.py --check
```

A healthy setup looks like this:

```
ldss3_qa 1.1.0 -- environment check

  [ok] python >= 3.11                     3.12.14
  [ok] PypeIt importable                  2.0.1.dev1241+gf3a1f1d27
  [ok] magellan_ldss3 present
  [ok] fork: joins the c1/c2 files        amp_files()
  [ok] fork: read noise by readout mode   ronoise_by_speed
  [ok] detector saturation sane           96,336 e-

  Ready.
```

Each line is checked for a reason. The two `fork:` lines matter most: stock
PypeIt *also* ships a `magellan_ldss3` module, but it cannot join the two
amplifier files and will silently reduce every exposure twice. This check tells
the two apart.

---

## Run it

```bash
python ldss3_qa.py /path/to/magellan_ldss3_A --tar
```

The argument is your **reduction** directory — the one holding `Calibrations/`
and `Science/` — not the directory of raw files.

| Option | Effect |
|---|---|
| `-o DIR` | write the bundle elsewhere (default: `<redux>/LDSS3_QA`) |
| `--tar` | also write `LDSS3_QA.tar.gz`, ready to send |
| `--no-pypeit-qa` | skip copying PypeIt's own QA folder, for a much smaller bundle |
| `--check` | check the environment and exit |
| `--version` | print the version |

It takes a minute or two, printing progress as it goes. The summary is printed
at the end as well as written into the bundle.

If part of the reduction is missing — no flats, say, or nothing extracted — the
affected checks are skipped with a note rather than failing, so it is safe to
run on a partial reduction.

---

## What you get

```
LDSS3_QA/
  SUMMARY.txt                     read this first
  <calib_key>/                    one folder per calibration set
    01_seam_continuity.png        flux match where the two amplifiers join
    02_slit_map.png               traced slits, amplifier seam marked
    03_wavecal_quality.png        dispersion, line count and RMS per slit
    03_wavecal_table.txt          the same, as numbers
    07_slit_saturation.png        which apertures are saturated alignment holes
  science/
    04_spec2d_overview.png        sky-subtracted 2D, slit and object traces
    05_wavelength_across_seam.png wavelength continuity across the join
    06_spec1d.png                 extracted spectra
    08_noise_model.png            is the variance model correct?
  pypeit_qa/                      PypeIt's own QA, copied unchanged
  *.pypeit *.calib *.par *.log    your inputs, so the run can be reproduced
```

---

## Reading the results

Open `SUMMARY.txt` first — it ends with a short checklist. In order of
importance:

### 1. Anything saying SUSPECT, OUTLIER or FAILED

A slit is flagged **SUSPECT** when its dispersion differs from the median of the
other slits by more than 5 %. Dispersion is the reliable discriminator here.
**Fit RMS is not** — a badly wrong solution can still show a small RMS if it
latched onto only a handful of lines.

Cross-check in `03_wavecal_table.txt`: a healthy slit has 50–70 identified arc
lines, a broken one has 10–30.

### 2. Apertures flagged as alignment holes

On a slitmask these saturate in the flats and arcs. Their arc lines merge into
flat-topped blends, so their wavelength solutions cannot be recovered. **Do not
extract science from them.** This is expected behaviour, not a failure of the
reduction.

### 3. The noise model

`sigma(chi)` is the scatter of (data − sky − object) × √ivar over the columns
where the sky was actually modelled. **It should be close to 1.**

| Value | Meaning |
|---|---|
| ≈ 1.0 | the uncertainties on your spectra can be believed |
| > 1.15 | the variance model understates the noise; error bars are too small |
| < 0.9 | the variance model overstates the noise |

If it is high, read noise is the usual cause. LDSS3's `ENOISE` header card is
unreliable — it reports the Slow-mode value whatever readout mode was used —
which is why the fork takes read noise from the published table instead. On a
correct reduction this check returns about 0.99; with the header value it
returned 1.30.

### 4. The amplifier seam

Two numbers, both about the join between the two halves:

- **Flux step in the combined flat** — typically 1–3 %, a small residual in the
  relative gain of the two amplifiers. After normalising it should fall below
  about 0.5 %. Measured on reference data: 2.82 % → 0.04 %, and 1.62 % → 0.09 %.
- **Wavelength jump across the seam** — should be far below 0.1 Å. Measured on
  reference data: 0.003 Å and 0.018 Å.

If either is much larger, the two halves were not joined correctly and the
reduction should not be trusted. Say so straight away.

Both checks need illuminated data on either side of the boundary, so they only
work when a slit actually lies across it. On a longslit that is almost always
true. On a slitmask the slits sit wherever the targets are, and if none happens
to fall on the boundary the checks report *not measurable* rather than inventing
a number. That is not a problem: the join is a property of the detector, not of
your mask, so it is either right for every dataset or wrong for every dataset —
this particular one simply cannot confirm it.

### 5. The 2D overview

In `04_spec2d_overview.png`, slit traces should sit on the illuminated slits and
object traces on the objects. Watch for slits with no trace, or traces drifting
off their slit along the dispersion direction.

---

## Troubleshooting

| Symptom | Cause and fix |
|---|---|
| `error: a required package is missing` | Wrong environment. `conda activate` the one you reduced in. |
| `--check`: `PypeIt importable ... not installed` | Same cause — PypeIt is not importable here. |
| `--check`: `you appear to have stock PypeIt, not the LDSS3 fork` | Install the fork; see [`LDSS3_REDUCTION_MANUAL.md`](LDSS3_REDUCTION_MANUAL.md) §2. |
| `--check`: `your fork predates the read-noise correction` | `cd` into your PypeIt clone, `git pull`, and re-reduce. Results from before that fix have error bars roughly 30 % too small. |
| `does not look like a PypeIt reduction directory` | Wrong folder. Point at the directory containing `Calibrations/` and `Science/` — not the raw data, and not its parent. |
| `no Slits_*.fits.gz found` | The reduction never got as far as tracing slits. Check the PypeIt `.log` for the real failure. |
| Everything is skipped | The reduction is incomplete. The summary names what was missing; the PypeIt log says where it stopped. |

---

## Sending results back

```bash
python ldss3_qa.py /path/to/magellan_ldss3_A --tar
```

Send `LDSS3_QA.tar.gz`. It carries the plots, `SUMMARY.txt`, and your `.pypeit`,
`.calib`, `.par` and `.log` files, so the run can be reproduced exactly rather
than guessed at.

Typically 50–150 MB with PypeIt's own QA included; add `--no-pypeit-qa` for a
few MB if that is awkward to send.

If something failed outright, send the PypeIt `.log` too, and say which command
you ran and what the last message on screen was.

---

## Commissioning report

`report/` holds a technical report on the LDSS3-C support in PypeIt, based on
the commissioning data used to test it. The built PDF is
`report/ldss3_pypeit_report.pdf`.

---

## Caveats

- Only 1×1 binning has been tested. The amplifier seam is assumed to lie at
  spatial pixel 1024, correct for unbinned data; for binned data that constant
  in `ldss3_qa.py` needs dividing by the spatial binning.
- The saturation check compares against the ADC ceiling (65535 ADU × gain, in
  electrons). The gain is read from your reduction, so it follows your readout
  mode. It is a flag, not a measurement.

---

## License

BSD 3-Clause. See [`LICENSE`](LICENSE).
