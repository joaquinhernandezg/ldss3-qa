# Reducing LDSS3 data with PypeIt — internal manual

Magellan/LDSS3-C, all three VPH grisms, longslit and multi-slit.
Fork: <https://github.com/joaquinhernandezg/PypeIt>, branch `ldss3_multigrating`.

This is an internal working document, not polished documentation. If something
here disagrees with what the pipeline actually does, the pipeline is right —
tell Joaquín.

---

## 1. What this fork changes

LDSS3-C reads one CCD through two amplifiers and writes **each amplifier to its
own file**:

```
ccd0042c1.fits     amplifier 1
ccd0042c2.fits     amplifier 2
```

This fork joins them inside the pipeline. **You do not merge them beforehand,
and you do not list both files in your `.pypeit` file.** `pypeit_setup` ingests
one file per exposure and finds the companion automatically.

After overscan subtraction and trimming the frame is 4096 × 2048, dispersion
along the first axis, with amplifier 1 occupying spatial pixels 0–1023 and
amplifier 2 occupying 1024–2047. That boundary is referred to below as the
**seam**.

It also adds VPH-Blue and VPH-Red wavelength solutions (upstream has VPH-All
only), frame typing suited to LDSS3 headers, a measured bad-column mask, and an
OH sky line list for sky-based wavelength calibration.

---

## 2. Install

The branch tracks PypeIt `develop`, which needs numpy ≥ 2.4 and astropy ≥ 7, so
give it its own environment.

```bash
conda create -y -n ldss3 python=3.12 pip
conda activate ldss3

# Install the Qt binding from conda first.  pip cannot build pyqt6 from source
# without a full Qt toolchain, and the install will fail if you skip this.
conda install -y -c conda-forge pyqt=6

git clone -b ldss3_multigrating https://github.com/joaquinhernandezg/PypeIt.git
cd PypeIt
pip install -e .
```

Check it took:

```bash
python -c "import pypeit; print(pypeit.__version__)"
pypeit_setup -h
```

The `-e` matters: it is an editable install, so `git pull` updates the code with
no reinstall.

### Keeping up to date

When Joaquín pushes changes:

```bash
cd PypeIt
git pull
```

That is all, unless the dependencies changed — if `git pull` brings in a new
`pyproject.toml`, re-run `pip install -e .`.

To pick up new upstream PypeIt work as well:

```bash
git remote add upstream https://github.com/pypeit/PypeIt.git   # once
git fetch upstream
git merge upstream/develop
```

If that conflicts, stop and ask rather than resolving it — the LDSS3 module is
the likely conflict and getting it wrong is silent.

---

## 3. Reducing

### Step 1 — look at what you have

```bash
cd /somewhere/with/space
pypeit_setup -r /path/to/raw -s magellan_ldss3
```

`-r` points at the directory holding the raw `c1`/`c2` pairs. This writes
`setup_files/magellan_ldss3.sorted`, listing every configuration it found and
how it typed each frame. **Read this before going further.**

A configuration is defined by grism, slitmask, binning and order-blocking
filter. You usually want the one with science frames plus arcs and flats.

### Step 2 — generate the input file

```bash
pypeit_setup -r /path/to/raw -s magellan_ldss3 -c A
```

Use the setup letter from the `.sorted` file. This makes
`magellan_ldss3_A/magellan_ldss3_A.pypeit`.

**Open it and check the `frametype` column.** This is the step people skip and
then wonder why the reduction is wrong. See §4.

### Step 3 — run

```bash
cd magellan_ldss3_A
run_pypeit magellan_ldss3_A.pypeit
```

Add `-o` to overwrite an earlier attempt. Expect tens of minutes to hours
depending on how many science frames you listed.

### Step 4 — QA

```bash
git clone https://github.com/joaquinhernandezg/ldss3-qa.git
python ldss3-qa/ldss3_qa.py magellan_ldss3_A --tar
```

Read `LDSS3_QA/SUMMARY.txt`, then send `LDSS3_QA.tar.gz` back. See §6.

---

## 4. Frame typing — the part to check

LDSS3's `EXPTYPE` header is not trustworthy on its own. Arcs and flats appear
under **both** `EXPTYPE='Flat'` and `EXPTYPE='Object'`, and bias frames report
whatever grism the wheel happened to be sitting on. So typing also reads
keywords out of the `OBJECT` card:

| `OBJECT` contains | Typed as | Real examples that match |
|---|---|---|
| `arc` `henear` `hene` `lamp` `comp` | `arc`, `tilt` | "arc ATM3a2_1", "LTT3218 arc", "GD71 VPH-RED OG590 Arc" |
| `flat` `quartz` `dome` `qh` `ff` | `pixelflat`, `illumflat`, `trace` | "flat ATM3a2_1", "EG274 - Qh", "…OG590 FF", "LTT3218 flatQh" |
| `zero` `bias` | `bias` | "Zero", "bias" |
| `dark` | `dark` | "dark" |
| `align` `thr` `field` `focus` `acq` | ignored | "align MASK …", "thr ATM3a2_1" |
| anything else | `science` or `standard` | "science ATM3a2_1", "EG274 - spec" |

`flat`, `quartz` and `dome` match anywhere in the name, so `flatQh` is caught.
The short ones must be whole words, so *Arcturus* is not an arc and an *offset*
frame is not a flat.

**Science vs standard is decided purely by exposure time** — longer or shorter
than 100 s by default. If your standards were long exposures, either edit the
`.pypeit` file by hand or set:

```ini
[calibrations]
    [[standardframe]]
        exprng = None, 300
[scienceframe]
    exprng = 300, None
```

Frames taken with the grism wheel open are acquisition or through-slit images
and are dropped at setup. That is intentional.

**If your `OBJECT` names follow none of these conventions**, the fix is simply
to edit the `frametype` column of the `.pypeit` file. That is a supported way to
work, not a hack.

---

## 5. Wavelength calibration

Arc lamps are the default for all three grisms:

| Grism | Typical coverage | Dispersion |
|---|---|---|
| VPH-Blue | 3660–6440 Å | 0.69 Å/px |
| VPH-Red | 5850–10590 Å | 1.18 Å/px |
| VPH-All | 3890–10470 Å | 1.93 Å/px |

Because arcs are taken at a different telescope pointing from the science
frames, residual spectral flexure is corrected against the sky by default
(`[flexure] spec_method = boxcar`). Leave that alone unless you switch to
sky-line calibration.

### Sky lines instead of arcs

For deep exposures, the sky OH forest can beat an arc-lamp solution, because it
is recorded through the same optical path at the same telescope orientation.
Most useful for VPH-Red. Add to the `.pypeit` file:

```ini
[calibrations]
    [[wavelengths]]
        lamps = OH_LDSS3_vac
        method = full_template
        reid_arxiv = magellan_ldss3_vph_red_sky.fits
    [[arcframe]]
        exprng = 0, None
    [[tiltframe]]
        exprng = 0, None
[flexure]
    spec_method = skip
```

and set the `frametype` of your science frames to `arc,tilt,science` so the sky
spectrum is used. Turning flexure off matters: a sky-derived solution is already
in the frame of the observation, so correcting it against the sky again
double-counts the shift.

The `OH_LDSS3_vac` list has 300 lines over 6000–10400 Å, in vacuum, from the
UVES sky atlas. It is useless blueward of 6000 Å where sky OH is sparse.

**This path is documented but has not been run end to end.** If you use it, say
so when you send the QA.

---

## 6. What to look at

Run `ldss3_qa.py` and read `SUMMARY.txt` top to bottom. In order of importance:

**1. Slits flagged SUSPECT.** A slit whose dispersion differs from the median of
the others by more than 5 % has a wrong solution. Trust dispersion, not fit RMS
— a badly wrong solution can still have a low RMS if it latched onto a handful
of lines. Check `03_wavecal_table.txt`: a good slit has 50–70 identified lines;
a broken one has 10–30.

**2. Slits flagged near saturation.** On a slitmask these are the **alignment
boxes**. They saturate in the flats and arcs, their arc lines merge into
flat-topped blends, and their wavelength solutions cannot be recovered. Do not
extract science from them. This is expected and is not something to fix.

**3. The seam.** Two numbers:

- Flux step in the combined flat: normally 1–3 %, a small residual in the
  relative gain of the two amplifiers. After normalising it should drop below
  about 0.5 %. Values measured in testing: 2.82 % → 0.04 %, and 1.62 % → 0.09 %.
- Wavelength jump across the seam: should be far below 0.1 Å. Measured:
  0.003 Å and 0.018 Å.

If either is much larger, the amplifier join is wrong and the reduction should
not be trusted. Say so immediately.

If no slit lies across the seam, both report "not measurable". That is normal
for some masks.

**4. The 2D overview** (`04_spec2d_overview.png`). Slit traces should sit on the
illuminated slits; object traces should sit on the objects. Look for slits with
no trace, or traces drifting off their slit.

---

## 7. Caveats — read before trusting anything

- **Only 1×1 binning has been tested.** Anything else logs a warning. The data
  and overscan sections are taken from the header and assumed already binned; if
  you have binned data, check the first frame carefully and tell Joaquín.

- **Alignment boxes are not excluded automatically.** PypeIt's
  `saturated_slits='mask'` needs more than half a slit's pixels saturated, and
  these do not reach that. They are traced and reduced as if they were science
  slits, and produce rubbish. Identify them from the QA and ignore them.

- **A slit whose edge lands on a bad column will be dropped.** About a dozen
  narrow bad-column groups are masked. If a slit you expect is missing, check
  its edges against the ranges in `bpm()` in
  `pypeit/spectrographs/magellan_ldss3.py`.

- **The order-blocking filter is part of the configuration.** If the filter was
  changed partway through a night without re-taking calibrations, the affected
  frames land in their own setup with no arcs or flats. That is deliberate —
  flat-fielding across a filter change would be worse — but it means those
  frames need a decision, not a rerun.

- **VPH-All is extrapolated blueward of ~3890 Å** (the first ~740 pixels). No
  arc line is detectable there and the grism has no throughput, so science
  frames carry no signal either. The wavelengths printed in that range are
  meaningless, not wrong data.

- **No flux calibration has been set up.** All spectra come out in counts.
  `pypeit_sensfunc` has not been run or validated for LDSS3.

- **VPH-Red and VPH-All are validated on longslit commissioning data only.**
  VPH-Blue is the only grism tested on a real slitmask.

---

## 7b. Detector numbers

Published LDSS3-C properties, and how the pipeline uses them.

| Property | Value |
|---|---|
| Full well | ~205,000 e- (10% non-linear) |
| Non-linearity (1%) | > 175,000 e- |
| ADC | 16-bit, so 65535 ADU maximum |
| Readout modes | Slow, Fast, Turbo |

Read noise, e-, by mode and amplifier:

| Mode | amp 1 (C1) | amp 2 (C2) |
|---|---|---|
| Slow | 4.4 | 5.3 |
| Fast | 7.0 | 7.2 |
| Turbo | ~10 | ~10 |

**Saturation is set by the ADC, not by the full well.** At a gain of
1.65 e-/ADU the most a pixel can record is 65535 x 1.65 = 108,133 e-, and at
1.47 it is 96,336 e-. Both are well below the 175,000 e- where the CCD starts
to depart from linear, so the response is linear right up to the point the
converter clips. A saturation level taken from the full well would sit above
anything the detector can produce and would never flag a single pixel. The
pipeline therefore sets it from the ADC ceiling.

**A caution on read noise.** Every frame we have seen carries `SPEED='Fast'`
in the header but `ENOISE` = 4.67 / 5.06 e-, which are the *Slow* values from
the table above; the Fast values are 7.0 / 7.2. The header does not appear to
track the readout mode. The pipeline uses the header values, because they are
per-frame and there is nothing better to use, but if the noise looks optimistic
in your reduction this is the first thing to check. Tell Joaquin if you see it.

Gain is read per-frame from `EGAIN` (1.65 and 1.47 e-/ADU in our data).

---

## 8. Useful extra settings

Cut spurious detections near slit edges on multi-slit data:

```ini
[reduce]
    [[findobj]]
        find_trim_edge = 100, 100
```

Inspect results interactively:

```bash
pypeit_show_2dspec Science/spec2d_*.fits        # 2D, needs ginga
pypeit_show_1dspec Science/spec1d_*.fits        # 1D
pypeit_chk_wavecalib Calibrations/WaveCalib_*.fits   # the per-slit table
pypeit_chk_edges Calibrations/Edges_*.fits.gz   # slit tracing
```

---

## 9. Sending results back

```bash
python ldss3-qa/ldss3_qa.py magellan_ldss3_A --tar
```

Send `LDSS3_QA.tar.gz`. It contains the plots, `SUMMARY.txt`, and your
`.pypeit`, `.calib`, `.par` and `.log` files, so the run can be reproduced
exactly. With PypeIt's own QA it is usually 50–150 MB; add `--no-pypeit-qa` for
a few MB if that is awkward to send.

If something failed outright, send the `.log` file too, and say which command
you ran and what the last thing on screen was.
