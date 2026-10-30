#!/usr/bin/env python3
"""Adds a side box to every pdf under plots/6_cuts/<cut>/ (called by plot.sh after tidy_plots.py).

The box lists the numbers of the plot's sample from 6_cuts/rank_FOM_theta.csv (written by Plot3D.cpp, every column of every cut folder is in it):
CE_recon_proton_in_cut/recon_proton_in_cut, CE_recon_proton_in_cut/CE_recon_before_cut, nonCE_recon_proton_in_cut/recon_proton_in_cut, FOM_theta, FOM_thetaBlock.
The page is widened to the right so the box never covers the plot. Plots of 0_KE_KECUT80_* and of 2_KECUT80 use the "cut and ZDC KE > 80%" row, the others the "cut" row.
A folder that holds two cuts (4_nProtons) gets no box on the plots that sit beside the cut folders, they belong to neither cut. pdfs already boxed get their box taken off and redone (unbox), so the box can change without rerunning the macros.
"""
import csv
import io
import os
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

from PIL import Image
from pypdf import PageObject, PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, DictionaryObject, NameObject, RectangleObject

COLS = ["CE_recon_proton_in_cut/recon_proton_in_cut", "CE_recon_proton_in_cut/CE_recon_before_cut",
        "nonCE_recon_proton_in_cut/recon_proton_in_cut", "FOM_theta", "FOM_thetaBlock"]
MARK = "/CutSideBox"
A4 = (595, 842) # pt, the page ROOT writes (portrait, /Rotate 90 for the landscape canvases)
WIDTH, FONT, LEAD = 270.0, 7.5, 11.0
BOX_TOP = 0.53 # box top as a fraction of the page height when the plot has no legend: where the legends of the 1D plots end (NDC y 0.44 of the 20 x 16 cm figure ROOT puts at the top of the A4 page)
LEGEND_W = (100, 200) # pt, width of the legend frame (0.245 of the 567 pt wide figure) used to tell it from the axes and the histogram lines
GAP = 12.0 # pt between the legend frame and the first line of the box
KE80 = "cut and ZDC KE > 80%"
FOM_CUTS = {"1_reconTheta": "1_hitZDC", "2_reconBlockTheta": "3_ResCuts/2mmCutResD"} # 5_fomCalc folder -> rank row cut of its sample


def load_rows(cuts_dir):
    rows = {}
    with open(os.path.join(cuts_dir, "rank_FOM_theta.csv"), newline="") as f:
        for r in csv.DictReader(f):
            cut = "1_hitZDC" if r["cut"].startswith("1_hitZDC") else r["cut"]
            rows[(cut, r["sample"])] = r
    return rows


def cut_of(rel_dir, rows):
    """The rank row key's cut for a pdf folder (relative to 6_cuts), or None when no single cut owns it."""
    parts = rel_dir.split("/")
    if parts[0] in ("1_hitZDC", "2_KECUT80"):
        return "1_hitZDC"
    if len(parts) >= 2 and (parts[0] + "/" + parts[1], "cut") in rows:
        return parts[0] + "/" + parts[1]
    sub = sorted({c for c, _ in rows if c.split("/")[0] == parts[0] and "/" in c})
    return sub[0] if len(sub) == 1 else None


def lines_for(row):
    """Text lines of the box, N_{...} is a subscript (see text_ops). The counts come first, then the ratios and the FOMs of COLS."""
    n = {"ce": int(float(row["CE_recon_proton_in_cut"])), "non": int(float(row["nonCE_recon_proton_in_cut"]))}
    out = [f"N_{{CE,recon,in cut}} = {n['ce']}", f"N_{{nonCE,recon,in cut}} = {n['non']}"]
    label = {COLS[0]: "N_{CE,recon,in cut} / N_{recon,in cut}", COLS[1]: "N_{CE,recon,in cut} / N_{CE,recon}", COLS[2]: "N_{nonCE,recon,in cut} / N_{recon,in cut}", COLS[3]: "FOM_{theta}", COLS[4]: "FOM_{thetaBlock}"}
    for c in COLS:
        out.append(f"{label.get(c, c)} = {float(row[c]):.4g}")
    return out


def esc(s):
    return s.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def text_ops(s, size):
    """Text operators of s with N_{sub} drawn as a smaller, lowered subscript (the font size and the rise are reset after it)."""
    ops = []
    for i, part in enumerate(s.split("_{")):
        if i == 0:
            ops.append(f"({esc(part)}) Tj")
            continue
        sub, _, rest = part.partition("}")
        ops.append(f"/FcsB {size * 0.8:.2f} Tf -{size * 0.25:.2f} Ts ({esc(sub)}) Tj /FcsB {size} Tf 0 Ts ({esc(rest)}) Tj")
    return " ".join(ops)


def legend_bottom(path):
    """(left edge, bottom edge) of the legend frame in pt from the left and the top of the displayed page, None without a legend. The page is rendered at 72 dpi (1 px = 1 pt, /Rotate honoured);
    a frame shows as two horizontal black runs of the same x extent and of the legend width (the stats, fit and CE boxes have it too). The legend is the frame with colored line samples
    inside. Without a legend the box goes under the lowest frame edge in the right half of the page (the stats or fit box); None when there is none either, and the box then goes at the right of the page."""
    png = subprocess.run(["gs", "-q", "-dBATCH", "-dNOPAUSE", "-dSAFER", "-sDEVICE=png16m", "-r72", "-dTextAlphaBits=1", "-dGraphicsAlphaBits=1", "-sOutputFile=-", path],
                         check=True, capture_output=True).stdout
    img = Image.open(io.BytesIO(png)).convert("RGB")
    w, h = img.size
    px = img.load()
    dark = lambda x, y: max(px[x, y]) < 100
    runs = {} # (x start, x end) -> [(row, x start, x end)]
    for y in range(h):
        x = 0
        while x < w:
            if not dark(x, y):
                x += 1
                continue
            x0 = x
            while x < w and dark(x, y):
                x += 1
            if LEGEND_W[0] <= x - x0 <= LEGEND_W[1]:
                runs.setdefault((round(x0 / 3), round(x / 3)), []).append((y, x0, x)) # rounded: the two edges differ by a pixel at most
    frames = [] # left, right, top, bottom
    for rows in runs.values():
        edges = [r for i, r in enumerate(rows) if i == 0 or r[0] - rows[i - 1][0] > 2] # one edge per line, a thick line is several rows
        for top, bottom in zip(edges[0::2], edges[1::2]): # frames sharing the x extent stack: edges 0 and 1 are one frame, 2 and 3 the next
            if bottom[0] - top[0] >= 15:
                frames.append((top[1], top[2], top[0], bottom[0]))
    white = lambda x, y: min(px[x, y]) > 230
    colored = lambda x, y: max(px[x, y]) - min(px[x, y]) > 80
    def is_legend(f): # a white inside with colored line samples, the dark cells of a 2D plot make no such frame
        pts = [(x, y) for y in range(f[2] + 2, f[3] - 1) for x in range(f[0] + 2, f[1] - 1)]
        return sum(white(x, y) for x, y in pts) > 0.6 * len(pts) and any(colored(x, y) for x, y in pts)
    legend = [f for f in frames if is_legend(f)]
    if not legend:
        right = [f for f in frames if f[0] > w / 2] # stats and fit boxes, the plot frame is wider than LEGEND_W; a lone run (a curve along the x axis) is no frame
        if not right: # no legend, stats or fit box: at the usual height, right of whatever is drawn there (the plot frame, or a color palette and its labels)
            top = int(h * (1 - BOX_TOP))
            used = [x for x in range(w) if any(dark(x, y) for y in range(max(0, top - 10), min(h, top + 110)))] # the box is about 100 pt high
            return (used[-1] + 12, top - GAP) if used else None
        key = (round(right[0][0] / 3), round(right[0][1] / 3))
        return right[0][0], max(r[0] for r in runs[key]) # the lowest edge of that column, a title line splits the stats box into two frames
    f = max(legend, key=lambda f: f[3])
    stacked = runs[(round(f[0] / 3), round(f[1] / 3))] # boxes stacked under the legend with its x extent (the CE lost box of FitCore.h DrawCELost, the residual formula box): the box goes under the lowest
    return f[0], max(r[0] for r in stacked)


def is_boxed(reader):
    """The side box is marked in the metadata, and also found in the page content (a viewer that saves the pdf, like Preview, can drop the metadata)."""
    if (reader.metadata or {}).get(MARK) == "1":
        return True
    return any(b"/FcsB" in page["/Contents"].get_object().get_data() for page in reader.pages if "/Contents" in page)


def unbox(path):
    """Takes the side box off a boxed pdf again, so the box can be redone with new numbers or a new place without rerunning the macros: the box is the last
    part of the page content (it starts at the q that clips to the widened page, before its first /FcsB text), and the page goes back to the A4 ROOT wrote
    (crop_pdfs.py may have cut it). Returns False when the pdf is not boxed."""
    reader = PdfReader(path)
    boxed = is_boxed(reader)
    if not boxed and all(list(p.mediabox) == [0, 0, A4[0], A4[1]] for p in reader.pages):
        return False
    writer = PdfWriter(clone_from=reader)
    for page in writer.pages:
        if not boxed: # cropped by crop_pdfs.py in an earlier run but never boxed (5_fomCalc): only back to the A4, stamp measures the page ROOT wrote
            page.mediabox = RectangleObject((0, 0, A4[0], A4[1]))
            page.cropbox = RectangleObject((0, 0, A4[0], A4[1]))
            continue
        if int(page.get("/Rotate", 0)) % 360 != 90:
            raise SystemExit(f"{path}: /Rotate not 90, cannot take the box off")
        stream = page["/Contents"].get_object()
        data = stream.get_data()
        text = data.find(b"/FcsB")
        start = data.rfind(b"\nq\n", 0, text) + 1 if text >= 0 else -1
        if start <= 0:
            raise SystemExit(f"{path}: no side box found in the page content")
        stream.set_data(data[:start])
        if "/Annots" in page:
            del page["/Annots"] # markup a viewer saved into the pdf goes with the old box
        page.mediabox = RectangleObject((0, 0, A4[0], A4[1]))
        page.cropbox = RectangleObject((0, 0, A4[0], A4[1]))
    writer.add_metadata({MARK: "0"})
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        writer.write(f)
    os.replace(tmp, path)
    return True


def stamp(path, lines, title):
    reader = PdfReader(path)
    if is_boxed(reader):
        return False
    writer = PdfWriter()
    legend = legend_bottom(path)
    for src in reader.pages:
        w, h = float(src.mediabox.width), float(src.mediabox.height)
        rot = int(src.get("/Rotate", 0)) % 360 # ROOT writes landscape canvases as a portrait page rotated 90
        if rot not in (0, 90):
            raise SystemExit(f"{path}: /Rotate {rot} not handled")
        dw, dh = (h, w) if rot else (w, h) # the page as displayed
        ov = PageObject.create_blank_page(width=w, height=h + WIDTH) if rot else PageObject.create_blank_page(width=w + WIDTH, height=h) # grown on the displayed right
        content = []
        for i, ln in enumerate([title] + lines):
            boxLeft, boxTop = (legend[0], dh - legend[1] - GAP) if legend is not None else (dw + 8, dh * BOX_TOP) # right under the legend frame, same left edge; without a legend at the right of the page
            dx, dy = boxLeft, boxTop - LEAD * (i + (0 if i == 0 else 0.6)) # displayed position
            m = f"0 1 -1 0 {w - dy:.1f} {dx:.1f}" if rot else f"1 0 0 1 {dx:.1f} {dy:.1f}"
            content.append(f"BT /FcsB {FONT} Tf {m} Tm {text_ops(ln, FONT)} ET")
        font = DictionaryObject({NameObject("/Type"): NameObject("/Font"), NameObject("/Subtype"): NameObject("/Type1"),
                                 NameObject("/BaseFont"): NameObject("/Helvetica")})
        ov[NameObject("/Resources")] = DictionaryObject({NameObject("/Font"): DictionaryObject({NameObject("/FcsB"): writer._add_object(font)})})
        stream = DecodedStreamObject()
        stream.set_data("\n".join(content).encode("latin-1"))
        ov[NameObject("/Contents")] = writer._add_object(stream)
        top, right = (src.mediabox.top + WIDTH, src.mediabox.right) if rot else (src.mediabox.top, src.mediabox.right + WIDTH)
        src.mediabox.upper_right = (right, top)
        src.cropbox.upper_right = (right, top)
        src.merge_page(ov)
        writer.add_page(src)
    writer.add_metadata({MARK: "1"})
    tmp = path + ".tmp"
    with open(tmp, "wb") as f:
        writer.write(f)
    os.replace(tmp, path)
    return True


def unbox_and_stamp(job):
    path, lines, title = job
    unbox(path) # always redone: a boxed pdf from an earlier run gets the current numbers and place
    return stamp(path, lines, title)


def main(plots):
    cuts_dir = os.path.join(plots, "6_cuts")
    rows = load_rows(cuts_dir)
    skipped = 0
    jobs = [] # (path, lines, title): collected first, then run in parallel, one worker per pdf
    for d, _, files in os.walk(cuts_dir):
        rel = os.path.relpath(d, cuts_dir)
        for f in sorted(files):
            if not f.endswith(".pdf"):
                continue
            cut = cut_of(rel, rows) if rel != "." else None
            ke80 = f.startswith("0_KE_KECUT80_") or rel.split(os.sep)[0] == "2_KECUT80"
            row = rows.get((cut, KE80 if ke80 else "cut"))
            if row is None:
                skipped += 1
                continue
            title = f"{cut if cut != '1_hitZDC' else '1_hitZDC (no cut)'}, {'ZDC KE > 80%' if ke80 else 'no KE cut'}"
            jobs.append((os.path.join(d, f), lines_for(row), title))
    for sub, cut in FOM_CUTS.items(): # 5_fomCalc samples (ke_gt80, PlotFomPolarization.cpp): recon theta has no cut, recon block theta the block residual cut (kBlockCut, QualityCut.h)
        for d, _, files in os.walk(os.path.join(plots, "5_fomCalc", sub)):
            for f in sorted(files):
                if f.endswith(".pdf"):
                    jobs.append((os.path.join(d, f), lines_for(rows[(cut, "cut")]), f"{cut if cut != '1_hitZDC' else '1_hitZDC (no cut)'}, no KE cut"))
    if len(jobs) > 1:
        with ProcessPoolExecutor(max_workers=os.cpu_count()) as ex:
            done = sum(ex.map(unbox_and_stamp, jobs, chunksize=4))
    else:
        done = sum(unbox_and_stamp(j) for j in jobs)
    print(f"cut side box: {done} pdfs boxed, {skipped} without a single owning cut")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "plots")
