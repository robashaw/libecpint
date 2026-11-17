#!/usr/bin/env python3
"""Compare ecpint_bench results (see README.md).

Reads CSV files written by `ecpint_bench --csv FILE --tag TAG`, and compares each tag with the
first (or --base). Several files with the same tag are repeated runs (rounds) of one build: each
measurement's median over them is used, and the spread between them counts towards the noise.
Run the builds alternately rather than one after another, so that drift (e.g. thermal throttling)
affects them alike.

For each measurement, prints the base median and, for each other tag, the ratio to it (below 1 is
faster) with the noise of the two: the larger of the spread of the per-run medians and the spread
of the repetitions within a run. Geometric means summarise the off-centre kernels, the kernels
with a shell on the ECP (controls for changes to the off-centre code), the systems and the
scaling. A `!` marks a checksum that differs from the base's: the two compute different values.

usage: compare.py FILE.csv [FILE.csv ...] [--base TAG] [--out FILE]
"""
import argparse
import csv
import math
import statistics
import sys


def checksums_agree(a, b):
    return abs(a - b) <= 1e-10 * max(abs(a), abs(b)) + 1e-14


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("files", nargs="+", help="CSV files from ecpint_bench --csv")
    p.add_argument("--base", help="tag to compare with (default: the first tag read)")
    p.add_argument("--out", help="also write the comparison to this CSV file")
    args = p.parse_args()

    # medians[key][tag] and spreads[key][tag]: one entry per run; key = (group, name)
    tags, order, medians, spreads, units, sums = [], [], {}, {}, {}, {}
    for path in args.files:
        with open(path) as f:
            for row in csv.DictReader(f):
                tag, key = row["tag"], (row["group"], row["name"])
                if tag not in tags:
                    tags.append(tag)
                if key not in medians:
                    order.append(key)
                    medians[key], spreads[key], sums[key] = {}, {}, {}
                medians[key].setdefault(tag, []).append(float(row["median"]))
                spreads[key].setdefault(tag, []).append(float(row["spread"]))
                units[key] = row["unit"]
                sums[key][tag] = float(row["checksum"])
    base = args.base or tags[0]
    if base not in tags:
        sys.exit("no results with tag %s" % base)
    others = [t for t in tags if t != base]
    runs = {t: max(len(medians[k].get(t, [])) for k in order) for t in tags}

    def summary(key, tag):
        values = medians[key][tag]
        m = statistics.median(values)
        between = (max(values) - min(values)) / m if m > 0 else 0.0
        return m, max(between, statistics.median(spreads[key][tag]))

    print("Runs per tag: " + ", ".join("%s %d" % (t, runs[t]) for t in tags))
    print("Ratio to %s (< 1 is faster), with the noise of the two (larger of the spread between"
          " runs and within runs)" % base)
    head = "  %-34s %12s %-7s" % ("", base, "")
    for t in others:
        head += " %18s" % t
    print(head)
    out_rows, logs, differ = [], {}, False
    for key in order:
        group, name = key
        if base not in medians[key]:
            continue
        bm, bs = summary(key, base)
        line = "  %-34s %12.4g %-7s" % (name, bm, units[key])
        for t in others:
            if t not in medians[key]:
                line += " %18s" % "-"
                continue
            m, s = summary(key, t)
            ratio = m / bm if bm > 0 else float("nan")
            agree = checksums_agree(sums[key][t], sums[key][base])
            differ |= not agree
            line += " %9.3f (%4.1f%%)%s" % (ratio, 100 * max(s, bs), "" if agree else " !")
            kind = group if group != "kernels" else (
                "kernels off" if name.startswith("off") else "kernels controls")
            if bm > 0 and m > 0:
                logs.setdefault((kind, t), []).append(math.log(ratio))
            out_rows.append([group, name, units[key], base, bm, t, m, ratio, max(s, bs),
                             sums[key][base], sums[key][t]])
        print(line)

    print("\nGeometric mean ratios")
    for kind in ["kernels off", "kernels controls", "systems", "scaling"]:
        cells = ["%9.3f (%s)" % (math.exp(sum(logs[(kind, t)]) / len(logs[(kind, t)])), t)
                 for t in others if logs.get((kind, t))]
        if cells:
            print("  %-34s %s" % (kind, "  ".join(cells)))
    if differ:
        print("\n! checksums differ from %s: the builds compute different values" % base)

    if args.out:
        with open(args.out, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["group", "name", "unit", "base", "base_median", "tag", "median", "ratio",
                        "noise", "base_checksum", "checksum"])
            w.writerows(out_rows)


if __name__ == "__main__":
    sys.exit(main())
