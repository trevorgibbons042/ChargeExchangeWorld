#!/usr/bin/env python3
"""Tidies plots/ after the plot macros ran (called last by plot.sh; the macros read and append to their csv files while they run).

1. Cut folders (ke_*, ke_cut_*) side by side in one area: their csv files with the same name and header become one csv in the area's
   csv/ folder, with a first column "cut" naming the folder (ke_gt80/1_apCE/stats.csv -> <area>/csv/1_apCE_stats.csv, cut = ke_gt80).
2. Any other folder holding more than one csv file gets them moved into its own csv/ subfolder (folders already named *csv stay).
3. A folder holding a single item (a plot or a folder) is removed and the item moves up into its parent, named <folder>_<item>
   (a csv/ folder with one file: the file keeps its name).
"""
import csv
import os
import re
import shutil
import sys

CUT_DIR = re.compile(r"^ke_(all|gt\d+|cut_\d+)(_res2mm)?$")


def csv_files(top):
    """(relative path without csv/ segments, absolute path) of every csv under top."""
    for d, _, files in os.walk(top):
        for f in files:
            if f.endswith(".csv"):
                rel = os.path.relpath(os.path.join(d, f), top)
                yield "/".join(p for p in rel.split(os.sep) if p != "csv"), os.path.join(d, f)


def merge_cut_folders(area):
    cuts = sorted(c for c in os.listdir(area) if CUT_DIR.match(c) and os.path.isdir(os.path.join(area, c)))
    if len(cuts) < 2:
        return
    groups = {}
    for cut in cuts:
        for rel, path in csv_files(os.path.join(area, cut)):
            groups.setdefault(rel, []).append((cut, path))
    for rel, members in groups.items():
        tables = []
        for cut, path in members:
            with open(path, newline="") as f:
                rows = list(csv.reader(f))
            if rows:
                tables.append((cut, rows))
        headers = {tuple(rows[0]) for _, rows in tables}
        if len(headers) != 1:
            continue  # different tables under one name: leave them where they are
        # numbers dropped from the path parts (1_theta/stats.csv -> theta_stats.csv): one number per folder, and the csv/ folder has no order
        out = os.path.join(area, "csv", "_".join(re.sub(r"^\d+_", "", p) for p in rel.split("/")))
        os.makedirs(os.path.dirname(out), exist_ok=True)
        with open(out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["cut"] + list(tables[0][1][0]))
            for cut, rows in tables:
                for r in rows[1:]:
                    w.writerow([cut] + r)
        for _, path in members:
            os.remove(path)


def remove_empty(top):
    for d, _, _ in sorted(os.walk(top), key=lambda x: -x[0].count(os.sep)):
        if d != top and not os.listdir(d):
            os.rmdir(d)


def gather_loose_csv(top):
    for d, _, files in os.walk(top):
        if os.path.basename(d).endswith("csv") or os.path.basename(d) == "6_cuts":  # the rank_*.csv files stay in 6_cuts
            continue
        loose = [f for f in files if f.endswith(".csv")]
        if len(loose) > 1:
            os.makedirs(os.path.join(d, "csv"), exist_ok=True)
            for f in loose:
                shutil.move(os.path.join(d, f), os.path.join(d, "csv", f))


def lift_single_items(top):
    changed = True
    while changed:
        changed = False
        for d, dirs, files in sorted(os.walk(top), key=lambda x: -x[0].count(os.sep)):
            if d == top or not os.path.isdir(d):
                continue
            items = os.listdir(d)
            if len(items) != 1 or os.path.basename(d) == "6_cuts" or "6_cuts" in d.split(os.sep) or CUT_DIR.match(os.path.basename(d)):  # everything under 6_cuts keeps its folders (1_hitZDC, 4_nProtons, ...)
                continue
            item = items[0]
            parent, name = os.path.split(d)
            target = os.path.join(parent, item if name == "csv" else name + "_" + item)  # a lone csv file needs no csv/ folder
            if os.path.exists(target):
                continue
            shutil.move(os.path.join(d, item), target)
            os.rmdir(d)
            changed = True


def lift_ke_into(parent):
    for name in ("ke_gt80", "fom_polarization"): # a fresh file replaces the old one of the same name up in the folder
        keDir = os.path.join(parent, name)
        if not os.path.isdir(keDir):
            continue
        for item in os.listdir(keDir):
            src, dst = os.path.join(keDir, item), os.path.join(parent, item)
            if os.path.isdir(src) and os.path.isdir(dst):
                shutil.copytree(src, dst, dirs_exist_ok=True)
                shutil.rmtree(src)
            elif os.path.isdir(src) or not os.path.exists(dst):
                shutil.move(src, dst)
            else:
                os.replace(src, dst)
        if not os.listdir(keDir):
            os.rmdir(keDir)


def lift_ke(top):
    """5_fomCalc (each theta set) and 3_ReconCalcVals have no KE cuts: the content of their ke_gt80 folder moves up."""
    area = os.path.join(top, "5_fomCalc")
    if os.path.isdir(area):
        for name in os.listdir(area):
            lift_ke_into(os.path.join(area, name))
    lift_ke_into(os.path.join(top, "3_ReconCalcVals"))


def main(top):
    lift_ke(top)
    for d, dirs, _ in os.walk(top):
        merge_cut_folders(d)
    remove_empty(top)
    gather_loose_csv(top)
    lift_single_items(top)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "plots")
