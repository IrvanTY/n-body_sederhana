import csv
import os
import statistics
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "src"))
from make_data import PLANETS, DWARFS, MOONS

TIMING_TOL = 1.02

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STRONG_CSV = os.path.join(ROOT, "results", "strong.csv")
WEAK_CSV = os.path.join(ROOT, "results", "weak.csv")
FIG_DIR = os.path.join(ROOT, "fig")
POS_INIT = os.path.join(ROOT, "data", "strong8000.csv")
POS_FINAL = os.path.join(ROOT, "results", "snap_final_n8000.csv")
POS_OUT_INIT = os.path.join(ROOT, "fig", "initial_n8000_xy.png")
POS_OUT_FINAL = os.path.join(ROOT, "fig", "final_n8000_xy.png")

CORE_NAMES = ["Sun"] + [body[0] for body in PLANETS] + [body[0] for body in DWARFS] + [body[0] for body in MOONS]
PLANET_NAMES = set(body[0] for body in PLANETS)
DWARF_NAMES = set(body[0] for body in DWARFS)


def load(path):
    with open(path) as handle:
        return list(csv.DictReader(handle))


def aggregate(rows):
    by_p = {}
    for row in rows:
        by_p.setdefault(int(row["P"]), []).append(row)
    stats = {}
    for procs, group in by_p.items():
        totals = [float(r["t_total"]) for r in group]
        comps = [float(r["t_comp"]) for r in group]
        comms = [float(r["t_comm"]) for r in group]
        spread = lambda values: statistics.stdev(values) if len(values) > 1 else 0.0
        sums = [float(r["t_sum"]) for r in group]
        stats[procs] = {
            "n": int(group[0]["N"]),
            "trials": len(group),
            "total": (statistics.mean(totals), spread(totals)),
            "comp": (statistics.mean(comps), spread(comps)),
            "comm": (statistics.mean(comms), spread(comms)),
            "accounted": statistics.mean(sums),
        }
    return stats


def check_timing(stats, label):
    for procs in sorted(stats):
        total, _ = stats[procs]["total"]
        accounted = stats[procs]["accounted"]
        if accounted > total * TIMING_TOL:
            print("%s P=%d inconsistent: max-rank t_sum %.6f > t_total %.6f"
                  % (label, procs, accounted, total), file=sys.stderr)
            return False
    return True


def karp_flatt(speedup, procs):
    if procs <= 1:
        return 0.0
    return (1.0 / speedup - 1.0 / procs) / (1.0 - 1.0 / procs)


def load_xy_init(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle))
    return [(float(r["x"]), float(r["y"])) for r in rows]


def load_xy_final(path):
    out = []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if line:
                tokens = line.split(",")
                out.append((float(tokens[0]), float(tokens[1])))
    return out


def make_xy_colors(n):
    cmap_parts = []
    for name in ("tab20", "tab20b", "tab20c"):
        try:
            cmap_parts.append(plt.get_cmap(name))
        except ValueError:
            pass
    cols = []
    for cmap in cmap_parts:
        cols.extend([cmap(i) for i in range(cmap.N)])
    return cols[:n]


def plot_xy(data, core, title, outpath, xlim, ylim, colors):
    fig, ax = plt.subplots(figsize=(12, 7))
    xa = [p[0] for p in data[core:]]
    ya = [p[1] for p in data[core:]]
    ax.scatter(xa, ya, s=2, c="#808080", alpha=0.35, linewidths=0, zorder=1, rasterized=True)
    for i in range(core):
        x, y = data[i]
        if CORE_NAMES[i] == "Sun":
            ax.scatter([x], [y], s=140, c=[colors[i]], marker="*",
                       edgecolors="black", linewidths=0.6, zorder=5)
        elif CORE_NAMES[i] in PLANET_NAMES:
            ax.scatter([x], [y], s=55, c=[colors[i]], marker="o",
                       edgecolors="black", linewidths=0.6, zorder=4)
        elif CORE_NAMES[i] in DWARF_NAMES:
            ax.scatter([x], [y], s=30, c=[colors[i]], marker="D",
                       edgecolors="black", linewidths=0.5, zorder=3)
        else:
            ax.scatter([x], [y], s=20, c=[colors[i]], marker="^",
                       edgecolors="black", linewidths=0.4, zorder=2)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x (AU)")
    ax.set_ylabel("y (AU)")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    handles = [Line2D([0], [0], marker="*", color="w", markerfacecolor=colors[0],
                       markeredgecolor="black", markersize=11, label=CORE_NAMES[0])]
    for i in range(1, core):
        marker = "o" if CORE_NAMES[i] in PLANET_NAMES else ("D" if CORE_NAMES[i] in DWARF_NAMES else "^")
        size = 8 if CORE_NAMES[i] in PLANET_NAMES else 6
        handles.append(Line2D([0], [0], marker=marker, color="w", markerfacecolor=colors[i],
                               markeredgecolor="black", markersize=size, label=f"{i}:{CORE_NAMES[i]}"))
    handles.append(Line2D([0], [0], marker=".", color="w", markerfacecolor="#808080",
                           markersize=8, label=f"Asteroid ({len(data) - core})"))
    fig.legend(handles=handles, loc="center right", bbox_to_anchor=(0.99, 0.5),
               fontsize=6, ncol=2, title="Bodies", markerscale=1.0)
    fig.subplots_adjust(right=0.78, left=0.07, top=0.90, bottom=0.12)
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print("wrote " + outpath)


def plot_positions():
    init = load_xy_init(POS_INIT)
    final = load_xy_final(POS_FINAL)
    assert len(init) == len(final), (len(init), len(final))
    assert len(init) >= len(CORE_NAMES)
    core = len(CORE_NAMES)
    xs = [p[0] for p in init + final]
    ys = [p[1] for p in init + final]
    pad = 0.02 * max(max(xs) - min(xs), max(ys) - min(ys))
    xlim = (min(xs) - pad, max(xs) + pad)
    ylim = (min(ys) - pad, max(ys) + pad)
    print(f"N={len(init)} xlim={xlim} ylim={ylim}")
    colors = make_xy_colors(core)
    plot_xy(init, core, "Posisi awal N=8000 strong8000 (proyeksi XY, AU)", POS_OUT_INIT, xlim, ylim, colors)
    plot_xy(final, core, "Posisi akhir N=8000 strong8000 (proyeksi XY, AU)", POS_OUT_FINAL, xlim, ylim, colors)


def main():
    os.makedirs(FIG_DIR, exist_ok=True)
    rows = load(STRONG_CSV)
    stats = aggregate(rows)
    if not check_timing(stats, "strong"):
        return 1
    ps = sorted(stats)
    if stats[ps[0]]["trials"] < 3:
        print("need at least 3 trials per point", file=sys.stderr)
        return 1
    t1, s1 = stats[1]["total"]
    means = [stats[p]["total"][0] for p in ps]
    stds = [stats[p]["total"][1] for p in ps]
    speeds = [t1 / m for m in means]
    yerr = [s * ((s1 / t1) ** 2 + (st / m) ** 2) ** 0.5 if m and t1 else 0.0
            for s, st, m in zip(speeds, stds, means)]
    effs = [s / p for s, p in zip(speeds, ps)]
    flats = [karp_flatt(s, p) for s, p in zip(speeds, ps)]
    comp_means = [stats[p]["comp"][0] for p in ps]
    comm_means = [stats[p]["comm"][0] for p in ps]

    plt.figure()
    plt.errorbar(ps, speeds, yerr=yerr,
                 marker="o", capsize=3, label="measured ratio of means, 1-sigma sample")
    plt.plot(ps, ps, linestyle="--", label="ideal")
    plt.xscale("log", base=2)
    plt.xlabel("P")
    plt.ylabel("speedup S = T1 / Tp")
    plt.legend()
    plt.savefig(os.path.join(FIG_DIR, "speedup.png"), dpi=150)
    plt.close()

    plt.figure()
    plt.plot(ps, effs, marker="o")
    plt.xscale("log", base=2)
    plt.xlabel("P")
    plt.ylabel("efficiency E = S / P")
    plt.savefig(os.path.join(FIG_DIR, "efficiency.png"), dpi=150)
    plt.close()

    plt.figure()
    plt.plot(ps, flats, marker="o")
    plt.xscale("log", base=2)
    plt.xlabel("P")
    plt.ylabel("Karp-Flatt serial fraction e")
    plt.savefig(os.path.join(FIG_DIR, "karp_flatt.png"), dpi=150)
    plt.close()

    plt.figure()
    labels = [str(p) for p in ps]
    plt.bar(labels, comp_means, label="t_comp max-rank mean")
    plt.bar(labels, comm_means, bottom=comp_means, label="t_comm max-rank mean")
    plt.scatter(labels, means, color="black", marker="x", label="t_total mean")
    plt.xlabel("P")
    plt.ylabel("seconds (mean of trials)")
    plt.legend()
    plt.savefig(os.path.join(FIG_DIR, "breakdown.png"), dpi=150)
    plt.close()

    print("P N trials t_total_mean t_total_std speedup efficiency karp_flatt t_comp_mean t_comm_mean t_sum_mean")
    for p, s, e, f, m, sd, cc, cm in zip(ps, speeds, effs, flats, means, stds, comp_means, comm_means):
        print("%d %d %d %.3f %.3f %.3f %.3f %.4f %.3f %.3f %.3f"
              % (p, stats[p]["n"], stats[p]["trials"], m, sd, s, e, f, cc, cm, stats[p]["accounted"]))

    if os.path.exists(WEAK_CSV):
        weak = aggregate(load(WEAK_CSV))
        if not check_timing(weak, "weak"):
            return 1
        wps = sorted(weak)
        plt.figure()
        plt.errorbar(wps,
                     [weak[p]["total"][0] for p in wps],
                     yerr=[weak[p]["total"][1] for p in wps],
                     marker="o", capsize=3)
        plt.xscale("log", base=2)
        plt.xlabel("P (N = 1000 x P)")
        plt.ylabel("wall time seconds (mean of trials)")
        plt.savefig(os.path.join(FIG_DIR, "weak.png"), dpi=150)
        plt.close()
        print("weak scaling written")
    plot_positions()
    print("plots written to " + FIG_DIR)
    return 0


if __name__ == "__main__":
    sys.exit(main())
