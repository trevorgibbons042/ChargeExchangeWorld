#!/usr/bin/env python3
"""Crops every pdf under plots/ to its content (called by plot.sh after cut_side_box.py).

ROOT always writes an A4 page with the figure at 20 cm width in its top left corner, the rest of the page is blank. What is painted is read from the page content (pdfgeom.py, no ghostscript render),
the bounding box of the non-white things is found and the page box is cut to it plus PAD.
Cropping an already cropped pdf finds the same content box again, so rerunning changes nothing.
"""
import os
import sys
from concurrent.futures import ProcessPoolExecutor

from pypdf import PdfReader, PdfWriter
from pypdf.generic import RectangleObject

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))) # FitCore.h runs these with python3 -I, which leaves out the script folder
import pdfgeom

PAD = 8.0 # pt kept around the content


def content_box(page):
    """(xmin, ymin, xmax, ymax) of everything painted on the page in pdf user space (what pdfgeom.page_prims reads from the content), None if the page is empty."""
    boxes = []
    for p in pdfgeom.page_prims(page):
        if p[0] == "seg":
            boxes.append((min(p[2][0], p[3][0]), min(p[2][1], p[3][1]), max(p[2][0], p[3][0]), max(p[2][1], p[3][1])))
        elif p[0] != "blank": # the white background of ROOT's canvas is no content
            boxes.append(p[2])
    if not boxes:
        return None
    return min(b[0] for b in boxes), min(b[1] for b in boxes), max(b[2] for b in boxes), max(b[3] for b in boxes)


def crop(path):
    reader = PdfReader(path)
    if len(reader.pages) != 1:
        return False
    page = reader.pages[0]
    box = content_box(page)
    if box is None:
        return False
    x0, y0, x1, y1 = pdfgeom.media(page)
    rect = (max(x0, box[0] - PAD), max(y0, box[1] - PAD), min(x1, box[2] + PAD), min(y1, box[3] + PAD)) # the padding is the same on any /Rotate
    if abs(rect[2] - rect[0] - (x1 - x0)) < 1 and abs(rect[3] - rect[1] - (y1 - y0)) < 1:
        return False # already as tight as it gets
    writer = PdfWriter(clone_from=reader)
    for p in writer.pages:
        p.mediabox = RectangleObject(rect)
        p.cropbox = RectangleObject(rect)
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        writer.write(f)
    os.replace(tmp, path)
    return True


def main(plots):
    paths = []
    for d, _, files in os.walk(plots):
        for f in sorted(files):
            if f.endswith(".pdf"):
                paths.append(os.path.join(d, f))
    # each pdf is cropped by exactly one worker (parsing the content of a pdf dominates the time)
    if len(paths) > 1:
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
            done = sum(ex.map(crop, paths, chunksize=4))
    else:
        done = sum(crop(p) for p in paths)
    print(f"crop pdfs: {done} of {len(paths)} pdfs cropped to their content")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "plots")
