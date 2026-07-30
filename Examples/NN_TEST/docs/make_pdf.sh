#!/bin/sh
set -eu

cd "$(dirname "$0")"
pdflatex -interaction=nonstopmode -halt-on-error RADCAL_NN.tex
pdflatex -interaction=nonstopmode -halt-on-error RADCAL_NN.tex
