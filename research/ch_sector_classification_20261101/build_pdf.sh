#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p qa/build qa/render qa/texmf-var qa/texmf-config qa/font-cache
export TEXMFVAR="$PWD/qa/texmf-var"
export TEXMFCONFIG="$PWD/qa/texmf-config"
export XDG_CACHE_HOME="$PWD/qa/font-cache"
export TEXMF="{/usr/share/texlive/texmf-dist,/usr/share/texmf}"
export TEXINPUTS=".:/usr/share/texlive/texmf-dist/tex//:/usr/share/texmf/tex//:"
export TEXFORMATS="$PWD/qa/build:"
if [[ ! -f qa/build/xelatex.fmt ]]; then
 (cd qa/build && xetex -ini -etex -jobname=xelatex -interaction=nonstopmode '\input xelatex.ini' > format.log)
fi
xelatex -interaction=nonstopmode -halt-on-error -output-directory=qa/build note.tex > qa/build/build.log
xelatex -interaction=nonstopmode -halt-on-error -output-directory=qa/build note.tex >> qa/build/build.log
cp qa/build/note.pdf note.pdf
pdftoppm -png -r 110 note.pdf qa/render/page
