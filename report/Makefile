# Rebuild the LDSS3-C commissioning report.
#
#   make                 tables, figures and PDF from the existing measurements
#   make measurements    re-run the detector and calibration measurements
#   make reductions      re-run the PypeIt reductions (hours; needs the raw data)
#
# LDSS3_RAW and LDSS3_REDUX point at the raw data and the reductions
# (defaults in scripts/style.py).

PYTHON ?= python
LDSS3_RAW ?= /panoramix1/estudiantes/jhernandez/LDSS3_REDUCTION/report_2026/raw
LDSS3_REDUX ?= /panoramix1/estudiantes/jhernandez/LDSS3_REDUCTION/report_2026/redux
export LDSS3_RAW LDSS3_REDUX

.PHONY: all pdf tables figures measurements reductions clean

all: pdf

reductions:
	bash scripts/run_reductions.sh $(LDSS3_RAW) $(dir $(LDSS3_REDUX))
	mkdir -p pypeit_files && cp $(LDSS3_REDUX)/*/*.pypeit pypeit_files/

measurements:
	$(PYTHON) scripts/measure_detector.py $(LDSS3_RAW) results
	$(PYTHON) scripts/measure_calibs.py

tables:
	$(PYTHON) scripts/make_tables.py

figures: tables
	$(PYTHON) scripts/make_figures.py

pdf: figures
	latexmk -pdf -interaction=nonstopmode -halt-on-error ldss3_pypeit_report.tex

clean:
	latexmk -C ldss3_pypeit_report.tex
