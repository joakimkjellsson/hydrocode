#!/bin/bash

TEX="/usr/texbin/pdflatex -halt-on-error"
RST2TEX=/sw/bin/rst2latex.py
RST2HTML=/sw/bin/rst2html.py

echo " --- Make HTML --- "
$RST2HTML hydrocode.rst > hydrocode.html

echo " --- Make Latex and PDF --- "
$RST2TEX hydrocode.rst > hydrocode.tex
$TEX hydrocode.tex > tex.out

