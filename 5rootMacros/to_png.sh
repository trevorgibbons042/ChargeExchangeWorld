#!/bin/bash
# Renders every pdf under the folder (default plots) to a png next to it with ghostscript and deletes the pdf, one gs per core.
# Called last by plot.sh: the pdfs must be final (side box and crop done), cut_side_box.py and crop_pdfs.py read pdfs only.
DIR=${1:-plots}
NPROC=$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)
find "$DIR" -name '*.pdf' -print0 | xargs -0 -P "$NPROC" -I{} sh -c \
    'gs -q -dSAFER -dBATCH -dNOPAUSE -sDEVICE=png16m -r150 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -o "${1%.pdf}.png" "$1" && rm "$1"' _ {}
echo "to_png: $(find "$DIR" -name '*.png' | wc -l | tr -d ' ') pngs, $(find "$DIR" -name '*.pdf' | wc -l | tr -d ' ') pdfs left"
