#!/usr/bin/env python3
"""Merge table.json + ident.out into a ranked markdown table."""
import json, os, sys, math, argparse
HERE = os.path.dirname(os.path.abspath(__file__))

def etastr(D, r):
    num = "".join(f"\\eta_{{{d}}}^{{{e}}}" for d, e in zip(D, r) if e > 0)
    den = "".join(f"\\eta_{{{d}}}^{{{-e}}}" for d, e in zip(D, r) if e < 0)
    if not num: num = "1"
    return f"{num}/{den}" if den else num

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--table", default="table.json")
    ap.add_argument("--ident", default="ident.out")
    ap.add_argument("--out", default="table.md")
    ap.add_argument("--top", type=int, default=10000)
    a = ap.parse_args()
    rows = json.load(open(os.path.join(HERE, a.table)))
    idt = {}
    for fn in (a.ident, "ident_multi.out"):
        pth = os.path.join(HERE, fn)
        if not os.path.exists(pth): continue
        for l in open(pth):
            i, _, v = l.partition("|")
            try: ii = int(i.strip())
            except ValueError: continue
            v = v.strip()
            if fn.endswith("multi.out"):
                v = v.split("|", 1)[-1].strip()
            if v.startswith("UNIDENT") and ii in idt: continue
            if ii in idt and not idt[ii].startswith("UNIDENT"): continue
            idt[ii] = v
    out = []
    out.append("| # | $N$ | $M,\\chi,w$ | $t$ | $\\deg t$ | rec | $\\lambda_1$ | $\\lambda_2$ | $\\lambda_{\\min}$ | $c$ | $k$ | $\\sigma_p$ | score | budget | $a_n$ | period |")
    out.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, d in enumerate(rows[:a.top]):
        per = idt.get(i, "")
        if d["known"]: per = (d["known"] + "; " + per).strip("; ")
        sl = ",".join(f"{p}:{v}" for p, v in sorted(d["slopes"].items()) if abs(v) > 0.05) or "-"
        out.append("| %d | %d | %d,%d,%d | $%s$ | %d | (%d,%d) | %.4f | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            i, d["N"], d["M"], d["chi"], d["w"], etastr(d["divs"], d["etaq_t"]), d["tdeg"],
            d["r"]+1, d["D"], d["lam1"],
            ("%.4f" % d["lam2"]) if d["lam2"] is not None else "--",
            ("%.4f" % d["lam_min"]) if d.get("lam_min") is not None else "--",
            ("%.6g" % d["c"]) if d["c"] is not None else "--",
            str(d["k"]), sl,
            ("%.3f" % d["score"]) if d["score"] is not None else "--",
            ("%.3f" % d["budget"]) if d["budget"] is not None else "--",
            ",".join(d["a"][:6]), per))
    open(os.path.join(HERE, a.out), "w").write("\n".join(out) + "\n")
    print(f"wrote {a.out}: {len(rows)} rows")

if __name__ == "__main__":
    main()
