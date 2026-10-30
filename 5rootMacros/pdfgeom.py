"""Reads what a ROOT pdf draws straight from its page content, so crop_pdfs.py and cut_side_box.py need no ghostscript render.

page_prims(page) walks the content stream (q Q cm re m l c h S s f F f* B B* b b* n W, BT ET Tf Tm Td Tj TJ, rg RG g G) and returns the painted things in page
user space (the unrotated pdf coordinates), in drawing order:
  ("seg",  i, (x0, y0), (x1, y1), rgb)        one stroked line piece
  ("fill", i, (xmin, ymin, xmax, ymax), rgb)  a filled path, as its bounding box
  ("text", i, (xmin, ymin, xmax, ymax), rgb)  a string, the box from the font size and the Helvetica widths
Clip rectangles (re W n, the canvas clip ROOT sets before each text) cut the boxes. A white fill is kept as ("blank", ...) (the white inside of a legend or stats box), a white stroke or text is dropped.
Callers that look for what is visible skip "blank". The same page the old ghostscript render showed gives the same boxes, up to the glyph widths of the text.
to_display(page, x, y) turns user space into the page as displayed (x from the left, y from the top), honouring /Rotate and the media box origin.
"""
import re

TOKEN = re.compile(rb"\((?:\\.|[^\\()])*\)|\[|\]|/[^\s/\[\]()<>]*|[-+]?(?:\d+\.?\d*|\.\d+)|[A-Za-z*'\"]+", re.S)
# Helvetica advance widths per 1000 em for the ASCII 32 to 126 (the afm values): Times, ROOT's other font, is a little narrower, so a box is never too small
WIDTHS = [278, 278, 355, 556, 556, 889, 667, 191, 333, 333, 389, 584, 278, 333, 278, 278] + [556] * 10 + [278, 278, 584, 584, 584, 556, 1015,
          667, 667, 722, 722, 667, 611, 778, 722, 278, 500, 667, 556, 833, 722, 778, 667, 778, 722, 667, 611, 722, 667, 944, 667, 667, 611,
          278, 278, 278, 469, 556, 333, 556, 556, 500, 556, 556, 278, 556, 556, 222, 222, 500, 222, 833, 556, 556, 556, 556, 333, 500, 278,
          556, 500, 722, 500, 500, 500, 334, 260, 334, 584]
ESCAPE = re.compile(rb"\\([0-7]{1,3}|.)", re.S)


def string_em(tok):
    """Width in em of a pdf string token like (abc\\(d), from the table above, 556 for anything outside ASCII."""
    body = ESCAPE.sub(lambda m: bytes([int(m.group(1), 8)]) if m.group(1)[:1].isdigit() else m.group(1), tok[1:-1])
    return sum(WIDTHS[b - 32] if 32 <= b <= 126 else 556 for b in body) / 1000.0


def _mul(m, n):
    """Matrix product m x n of two pdf matrices (a b c d e f), m applied first."""
    return (m[0] * n[0] + m[1] * n[2], m[0] * n[1] + m[1] * n[3],
            m[2] * n[0] + m[3] * n[2], m[2] * n[1] + m[3] * n[3],
            m[4] * n[0] + m[5] * n[2] + n[4], m[4] * n[1] + m[5] * n[3] + n[5])


def _apply(m, x, y):
    return (m[0] * x + m[2] * y + m[4], m[1] * x + m[3] * y + m[5])


def _is_white(c):
    return min(c) > 0.98


def _box(pts):
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return (min(xs), min(ys), max(xs), max(ys))


def _clip(b, c):
    if c is None:
        return b
    r = (max(b[0], c[0]), max(b[1], c[1]), min(b[2], c[2]), min(b[3], c[3]))
    return r if r[0] <= r[2] and r[1] <= r[3] else None


def content_data(page):
    c = page["/Contents"].get_object()
    if hasattr(c, "get_data"):
        return c.get_data()
    return b"\n".join(x.get_object().get_data() for x in c)


def page_prims(page):
    ctm, fill, stroke = (1, 0, 0, 1, 0, 0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0)
    clip = None
    stack = []
    path = [] # subpaths, each a list of user space points; a re is stored as its four corners closed
    pending_clip = False
    tm = tlm = (1, 0, 0, 1, 0, 0)
    fsize = 0.0
    prims = []
    args = []
    n = 0

    def paint(do_fill, do_stroke):
        nonlocal n
        n += 1
        if do_fill and path:
            b = _clip(_box([p for sp in path for p in sp]), clip)
            if b:
                prims.append(("blank" if _is_white(fill) else "fill", n, b, fill))
        if do_stroke and not _is_white(stroke):
            for sp in path:
                for a, b in zip(sp, sp[1:]):
                    if clip is None or _clip(_box([a, b]), clip):
                        prims.append(("seg", n, a, b, stroke))

    def show(em):
        nonlocal n, tm
        n += 1
        if _is_white(fill) or not em or fsize <= 0:
            return
        advance = _mul((1, 0, 0, 1, em * fsize, 0), tm) # the next string of the line starts where this one ends (the subscripts of cut_side_box.py)
        m = _mul(tm, ctm)
        w, h = em * fsize, fsize
        b = _clip(_box([_apply(m, 0, -0.25 * h), _apply(m, w, -0.25 * h), _apply(m, w, 0.85 * h), _apply(m, 0, 0.85 * h)]), clip)
        if b:
            prims.append(("text", n, b, fill))
        tm = advance

    for tok in TOKEN.findall(content_data(page)):
        c = tok[:1]
        if c in b"([]/-+.0123456789" and tok not in (b"'", b'"'):
            args.append(tok)
            continue
        op = tok.decode("latin1")
        try:
            nums = [float(a) for a in args if a[:1] in b"-+.0123456789"]
        except ValueError:
            nums = []
        if op == "q":
            stack.append((ctm, fill, stroke, clip))
        elif op == "Q":
            if stack:
                ctm, fill, stroke, clip = stack.pop()
        elif op == "cm" and len(nums) == 6:
            ctm = _mul(tuple(nums), ctm)
        elif op == "rg" and len(nums) == 3:
            fill = tuple(nums)
        elif op == "RG" and len(nums) == 3:
            stroke = tuple(nums)
        elif op == "g" and len(nums) == 1:
            fill = (nums[0],) * 3
        elif op == "G" and len(nums) == 1:
            stroke = (nums[0],) * 3
        elif op == "re" and len(nums) == 4:
            x, y, w, h = nums
            path.append([_apply(ctm, x, y), _apply(ctm, x + w, y), _apply(ctm, x + w, y + h), _apply(ctm, x, y + h), _apply(ctm, x, y)])
        elif op == "m" and len(nums) == 2:
            path.append([_apply(ctm, *nums)])
        elif op in ("l", "c", "v", "y") and len(nums) >= 2 and path:
            path[-1].append(_apply(ctm, nums[-2], nums[-1]))
        elif op == "h" and path and len(path[-1]) > 1:
            path[-1].append(path[-1][0])
        elif op in ("S", "s"):
            paint(False, True)
            path = []
        elif op in ("f", "F", "f*"):
            paint(True, False)
            path = []
        elif op in ("B", "B*", "b", "b*"):
            paint(True, True)
            path = []
        elif op == "n":
            if pending_clip and path:
                clip = _clip(_box([p for sp in path for p in sp]), clip) or (0, 0, 0, 0)
            pending_clip = False
            path = []
        elif op in ("W", "W*"):
            pending_clip = True
        elif op == "BT":
            tm = tlm = (1, 0, 0, 1, 0, 0)
        elif op == "Tf" and len(nums) == 1:
            fsize = nums[0]
        elif op == "Tm" and len(nums) == 6:
            tm = tlm = tuple(nums)
        elif op == "Td" and len(nums) == 2:
            tm = tlm = _mul((1, 0, 0, 1, nums[0], nums[1]), tlm)
        elif op in ("Tj", "'", '"'):
            show(sum(string_em(a) for a in args if a[:1] == b"("))
        elif op == "TJ":
            show(sum(string_em(a) for a in args if a[:1] == b"("))
        args = []
    return prims


def media(page):
    b = page.mediabox
    return float(b.left), float(b.bottom), float(b.right), float(b.top)


def to_display(page, x, y):
    """User space to the page as displayed: x from the left, y from the top."""
    x0, y0, x1, y1 = media(page)
    rot = int(page.get("/Rotate", 0)) % 360
    if rot == 0:
        return x - x0, y1 - y
    if rot == 90:
        return y - y0, x - x0
    raise SystemExit(f"/Rotate {rot} not handled")


def display_size(page):
    x0, y0, x1, y1 = media(page)
    return (y1 - y0, x1 - x0) if int(page.get("/Rotate", 0)) % 360 == 90 else (x1 - x0, y1 - y0)


def prim_box(page, p):
    """Bounding box (left, top, right, bottom) of a prim as displayed."""
    pts = [p[2], p[3]] if p[0] == "seg" else [(p[2][0], p[2][1]), (p[2][2], p[2][3])]
    d = [to_display(page, *q) for q in pts]
    return (min(q[0] for q in d), min(q[1] for q in d), max(q[0] for q in d), max(q[1] for q in d))
