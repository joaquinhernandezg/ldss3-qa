#!/usr/bin/env bash
# Reproduce the reductions analysed in the commissioning report.
#
#   RAW   directory holding one sub-directory of *c1.fits/*c2.fits per night
#   WORK  directory for the reductions (needs ~20 GB)
#
# Usage: bash run_reductions.sh RAW WORK
set -euo pipefail
RAW=$(realpath "$1"); WORK=$(realpath -m "$2")
HERE=$(dirname "$(realpath "$0")")
mkdir -p "$WORK"; cd "$WORK"

# One pass over all raw data: writes setup_files/ and one .pypeit per setup
pypeit_setup -s magellan_ldss3 -r "$RAW"/*/ -v 1
pypeit_setup -s magellan_ldss3 -r "$RAW"/*/ -c A,B,C,D,F,G -v 1

# Curate the five datasets used in the report
python "$HERE/prepare_pypeit_files.py" "$WORK" "$WORK/redux"

# Reduce.  The VPH-Blue longslit data contain no science exposure, so only the
# calibrations are processed.
for ds in vph_all_longslit vph_red_longslit vph_blue_mos vph_red_mos; do
    (cd "$WORK/redux/$ds" && run_pypeit "$ds.pypeit" -o > run.log 2>&1) &
done
(cd "$WORK/redux/vph_blue_longslit" && run_pypeit vph_blue_longslit.pypeit -o -c > run.log 2>&1) &
wait
