#!/usr/bin/env python3
"""Crops every pdf under plots/ to its content (called by plot.sh after cut_side_box.py).

ROOT always writes an A4 page with the figure at 20 cm width in its top left corner, the rest of the page is blank. The page is rendered with
ghostscript at 72 dpi (1 px = 1 pt, /Rotate honoured), the non-white bounding box is found with PIL and the page box is cut to it plus PAD.
Cropping an already cropped pdf finds the same content box again, so rerunning changes nothing.
"""
import io
import os
import subprocess
import sys

from PIL import Image, ImageChops
from pypdf import PdfReader, PdfWriter
from pypdf.generic import RectangleObject

PAD = 8.0 # pt kept around the content


def content_box(path):
    """(left, top, right, bottom) in px = pt of the displayed page, None if the page is empty."""
    png = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-dSAFER", "-sDEVICE=png16m", "-r72", "-dTextAlphaBits=1", "-dGraphicsAlphaBits=1", "-sOutputFile=-", path],
                         check=True, capture_output=True).stdout
    img = Image.open(io.BytesIO(png)).convert("RGB")
    return img.size, ImageChops.difference(img, Image.new("RGB", img.size, (255, 255, 255))).getbbox()


def crop(path):
    (dw, dh), box = content_box(path)
    if box is None:
        return False
    left, top, right, bottom = max(0, box[0] - PAD), max(0, box[1] - PAD), min(dw, box[2] + PAD), min(dh, box[3] + PAD)
    reader = PdfReader(path)
    if len(reader.pages) != 1:
        return False
    page = reader.pages[0]
    w, h = float(page.mediabox.width), float(page.mediabox.height) # unrotated page
    rot = int(page.get("/Rotate", 0)) % 360
    if rot == 0: # displayed x = x, displayed y from the top = h - y
        rect = (left, h - bottom, right, h - top)
    elif rot == 90: # displayed (X, Y from the bottom) = (y, w - x)
        rect = (w - (dh - top), left, w - (dh - bottom), right)
    else:
        raise SystemExit(f"{path}: /Rotate {rot} not handled")
    if abs(rect[2] - rect[0] - float(page.mediabox.width)) < 1 and abs(rect[3] - rect[1] - float(page.mediabox.height)) < 1:
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
    done = total = 0
    for d, _, files in os.walk(plots):
        for f in sorted(files):
            if f.endswith(".pdf"):
                total += 1
                done += crop(os.path.join(d, f))
    print(f"crop pdfs: {done} of {total} pdfs cropped to their content")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "plots")
