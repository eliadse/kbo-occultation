# examples/make_plot_report.py
"""
Collect the analysis plots into one dated PDF report (for the supervisor log).

Sections are one per source folder; the run parameters (night, cutoff, detrend
and DC windows, bins, time range) are parsed from the filenames and printed
under each figure. Each run writes a new timestamped file, so older reports
are kept for comparison.

Usage:
    python examples/make_plot_report.py [--dirs DIR ...] [--dpi 150] [--out PATH]
"""
import argparse
import datetime as dt
import os
import re
import subprocess
from dataclasses import dataclass, field

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image

OUT_DIR = os.path.join("outputs", "reports")
PAGE = (11.69, 8.27)  # A4 landscape, inches
PORTRAIT = PAGE[::-1]  # for tall figures (stacked lightcurves)
IMG_BOX = [0.03, 0.17, 0.94, 0.74]  # figure-page image axes (left, bottom, width, height)


def _num(s):
    return float(s) if s is not None else -1.0


def _dc_combined(m):
    return {"Night": m["night"], "Cutoff": f"{m['cut']} Hz", "Bins": f"{m['bins']} samples"}


def _detrended(m):
    return {"Night": m["night"], "Detrend window": f"{m['win']} s", "DC window": f"{m['dc']} s",
            "Cutoff": m["cut"] or "none", "Bins": m["bins"],
            "Time range": "good-time window" if m["range"] == "goodtime" else "whole run"}


def _nb_combined(m):
    return {"Night": m["night"], "Bins": f"{m['bins']} samples", "Time range": m["range"]}


@dataclass
class Group:
    """One subsection: files matching `pattern` inside a section folder."""
    title: str
    pattern: re.Pattern = None             # None: match anything
    params: callable = None                # match -> {label: value}
    sort_key: callable = None              # match -> sortable
    by_night: bool = False                 # split into one subsection per night


@dataclass
class Section:
    title: str
    folder: str
    groups: list = field(default_factory=list)


SECTIONS = [
    Section("DC combined - frequency-cut scan", "outputs/dc_combine_cut_frequencies", [
        Group("Night", re.compile(r"combined_dc_(?P<night>\d{4}-\d\d-\d\d)_cut(?P<cut>[\d.]+)Hz_(?P<bins>\d+)sample_bins"),
              _dc_combined, lambda m: (_num(m["cut"]), int(m["bins"])), by_night=True),
    ]),
    Section("Detrended lightcurves", "outputs/detrended_lightcurves", [
        Group("Night", re.compile(r"detrended_(?P<night>\d{4}-\d\d-\d\d)_win(?P<win>[\d.]+)s_dc(?P<dc>[\d.]+)s"
                                  r"(?:_cut(?P<cut>[^_]+))?_(?P<bins>\d+)bins_(?P<range>goodtime|wholerun)"),
              _detrended, lambda m: (_num(m["win"]), _num(m["dc"]), m["cut"] or "", m["range"]), by_night=True),
    ]),
    Section("Notebook figures", "notebooks", [
        Group("Combined DC lightcurves and histograms",
              re.compile(r"Combined_dc_\w+?_(?P<night>\d{4}-\d\d-\d\d)_(?P<bins>\d+)sample_bins_(?P<range>.+)"),
              _nb_combined, lambda m: (int(m["bins"]), m["range"]), by_night=True),
        Group("Theoretical SNR model"),
    ]),
]


@dataclass
class Entry:
    path: str
    section: str
    subsection: str
    params: dict


def section_for(folder):
    norm = os.path.normpath(folder)
    for s in SECTIONS:
        if os.path.normpath(s.folder) == norm:
            return s
    return Section(os.path.basename(norm) or norm, folder, [Group("")])


def collect(sections):
    """Returns [(section, [(subsection, [Entry])])] in report order."""
    out = []
    for s in sections:
        files = sorted(f for f in os.listdir(s.folder) if f.lower().endswith(".png"))
        claimed, subs = set(), []
        for g in s.groups:
            hits = []
            for f in files:
                if f in claimed:
                    continue
                m = g.pattern.match(os.path.splitext(f)[0]) if g.pattern else None
                if g.pattern and not m:
                    continue
                claimed.add(f)
                hits.append((f, m))
            if g.sort_key:
                hits.sort(key=lambda h: g.sort_key(h[1]))
            buckets = {}
            for f, m in hits:
                sub = f"{g.title} {m['night']}" if g.by_night else g.title
                params = g.params(m) if g.params else {}
                buckets.setdefault(sub, []).append(Entry(os.path.join(s.folder, f), s.title, sub, params))
            subs.extend(sorted(buckets.items()) if g.by_night else buckets.items())
        # anything no group matched still goes in, filename-only caption
        rest = [f for f in files if f not in claimed]
        if rest:
            subs.append(("Other", [Entry(os.path.join(s.folder, f), s.title, "Other", {}) for f in rest]))
        out.append((s, subs))
    return out


def git_info():
    def run(*cmd):
        try:
            return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()
        except (OSError, subprocess.CalledProcessError):
            return ""
    commit = run("git", "rev-parse", "--short", "HEAD") or "unknown"
    branch = run("git", "rev-parse", "--abbrev-ref", "HEAD") or "unknown"
    dirty = bool(run("git", "status", "--porcelain"))
    return commit, branch, dirty


def footer(fig, page, stamp):
    fig.text(0.03, 0.02, f"Generated {stamp}", fontsize=8, color="0.4")
    fig.text(0.97, 0.02, f"Page {page}", fontsize=8, color="0.4", ha="right")


def title_page(pdf, stamp, layout):
    commit, branch, dirty = git_info()
    fig = plt.figure(figsize=PAGE)
    fig.text(0.5, 0.72, "KBO occultation - analysis plots", fontsize=26, ha="center", weight="bold")
    fig.text(0.5, 0.64, f"Created {stamp}", fontsize=14, ha="center")
    fig.text(0.5, 0.59, f"git {branch} @ {commit}" + ("  (uncommitted changes)" if dirty else ""),
             fontsize=11, ha="center", color="0.3", family="monospace")
    y = 0.48
    fig.text(0.25, y, "Source folders", fontsize=13, weight="bold")
    for s, subs in layout:
        y -= 0.045
        n = sum(len(e) for _, e in subs)
        fig.text(0.27, y, f"{s.folder}/", fontsize=11, family="monospace")
        fig.text(0.75, y, f"{n} plots", fontsize=11, ha="right")
    footer(fig, 1, stamp)
    pdf.savefig(fig)
    plt.close(fig)


def toc_page(pdf, stamp, toc):
    fig = plt.figure(figsize=PAGE)
    fig.text(0.08, 0.90, "Contents", fontsize=20, weight="bold")
    lines = len(toc)
    size = 11 if lines <= 30 else max(6, 11 * 30 / lines)
    step = min(0.035, 0.80 / max(lines, 1))
    y = 0.84
    for level, text, page in toc:
        x = 0.08 if level == 0 else 0.12
        fig.text(x, y, text, fontsize=size + (1 if level == 0 else 0), weight="bold" if level == 0 else "normal")
        fig.text(0.92, y, str(page), fontsize=size, ha="right")
        y -= step
    footer(fig, 2, stamp)
    pdf.savefig(fig)
    plt.close(fig)


def divider_page(pdf, stamp, page, number, section, subs):
    fig = plt.figure(figsize=PAGE)
    fig.text(0.5, 0.58, f"{number}. {section.title}", fontsize=24, ha="center", weight="bold")
    fig.text(0.5, 0.51, f"{section.folder}/", fontsize=12, ha="center", family="monospace", color="0.3")
    y = 0.42
    for sub, entries in subs:
        fig.text(0.5, y, f"{sub or 'Plots'}  ({len(entries)})", fontsize=12, ha="center")
        y -= 0.035
    footer(fig, page, stamp)
    pdf.savefig(fig)
    plt.close(fig)


def load_image(path, dpi):
    """Image downsampled to the resolution it is shown at, and the page size to show it on."""
    im = Image.open(path)
    im = im.convert("RGBA") if im.mode in ("P", "LA") else im
    size = PORTRAIT if im.height > im.width else PAGE
    max_w, max_h = int(size[0] * IMG_BOX[2] * dpi), int(size[1] * IMG_BOX[3] * dpi)
    scale = min(max_w / im.width, max_h / im.height, 1.0)
    if scale < 1.0:
        im = im.resize((round(im.width * scale), round(im.height * scale)), Image.LANCZOS)
    return np.asarray(im), size


def figure_page(pdf, stamp, page, number, entry, dpi):
    img, size = load_image(entry.path, dpi)
    fig = plt.figure(figsize=size)
    head = f"{number}. {entry.section}" + (f"  -  {entry.subsection}" if entry.subsection else "")
    fig.text(0.03, 0.955, head, fontsize=11, weight="bold")
    ax = fig.add_axes(IMG_BOX)
    ax.imshow(img, interpolation="none")
    ax.set_axis_off()

    name = os.path.basename(entry.path)
    caption = name if entry.params or entry.subsection != "Theoretical SNR model" \
        else os.path.splitext(name)[0].replace("_", " ")
    fig.text(0.03, 0.14, caption, fontsize=10, family="monospace")
    mtime = dt.datetime.fromtimestamp(os.path.getmtime(entry.path)).strftime("%Y-%m-%d %H:%M")
    fig.text(0.03, 0.122, f"file modified {mtime}", fontsize=8, color="0.4")
    if entry.params:
        tab = fig.add_axes([0.03, 0.045, 0.94, 0.07])
        tab.set_axis_off()
        t = tab.table(cellText=[list(entry.params.values())], colLabels=list(entry.params),
                      cellLoc="center", loc="center")
        t.auto_set_font_size(False)
        t.set_fontsize(9)
        for (r, _), cell in t.get_celld().items():
            if r == 0:
                cell.set_facecolor("0.9")
                cell.set_text_props(weight="bold")
    footer(fig, page, stamp)
    pdf.savefig(fig, dpi=dpi)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dirs", nargs="+", default=None, help="source folders (default: the three standard ones)")
    ap.add_argument("--dpi", type=int, default=150, help="resolution of the embedded images")
    ap.add_argument("--out", default=None, help="output PDF (default: outputs/reports/plot_report_<timestamp>.pdf)")
    args = ap.parse_args()

    now = dt.datetime.now()
    stamp = now.strftime("%Y-%m-%d %H:%M")
    out = args.out or os.path.join(OUT_DIR, f"plot_report_{now:%Y-%m-%d_%H%M}.pdf")
    os.makedirs(os.path.dirname(out) or ".", exist_ok=True)

    sections = [section_for(d) for d in args.dirs] if args.dirs else SECTIONS
    layout = collect(sections)

    # first pass: page numbers for the contents page
    toc, page = [], 3
    for i, (s, subs) in enumerate(layout, 1):
        toc.append((0, f"{i}. {s.title}", page))
        page += 1
        for sub, entries in subs:
            if sub:
                toc.append((1, sub, page))
            page += len(entries)

    meta = {"Title": "KBO occultation - analysis plots", "Author": "Elia do Souto", "CreationDate": now}
    with PdfPages(out, metadata=meta) as pdf:
        title_page(pdf, stamp, layout)
        toc_page(pdf, stamp, toc)
        page = 3
        for i, (s, subs) in enumerate(layout, 1):
            divider_page(pdf, stamp, page, i, s, subs)
            page += 1
            for _, entries in subs:
                for e in entries:
                    figure_page(pdf, stamp, page, i, e, args.dpi)
                    page += 1
    print(f"{out}: {page - 1} pages, {os.path.getsize(out) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
