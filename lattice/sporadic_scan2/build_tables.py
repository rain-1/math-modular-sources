#!/usr/bin/env python3
"""Build table.md (all rows) and top_table.md (top by budget, k>=1) from
table.json + ident.out, re-keying identifications by the canonical a-sequence."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from fmt_period import fmt

rows = json.load(open(os.path.join(HERE, "table.json")))
prevp = os.path.join(HERE, "table_prev.json")
prev = json.load(open(prevp)) if os.path.exists(prevp) else rows
idt = {}
ip = os.path.join(HERE, "ident.out")
if os.path.exists(ip):
    for l in open(ip):
        i, _, v = l.partition("|")
        try: ii = int(i.strip())
        except ValueError: continue
        if ii < len(prev): idt[tuple(prev[ii]["canon"])] = fmt(v)

def etastr(D, r):
    n = "".join(f"\\eta_{{{d}}}^{{{e}}}" for d, e in zip(D, r) if e > 0)
    dd = "".join(f"\\eta_{{{d}}}^{{{-e}}}" for d, e in zip(D, r) if e < 0)
    return (n or "1") + (f"/{dd}" if dd else "")

HDR = ("| # | $N$ | $M,\\chi,w$ | $t$ | $\\deg t$ | rec | $\\lambda_1$ | $\\lambda_2$ | "
       "$\\lambda_{\\min}$ | $c$ | $k$ | $\\sigma_p$ | score | budget | $a_n$ | period |")
SEP = "|" + "---|"*16

def line(i, d):
    per = idt.get(tuple(d["canon"]), "")
    if d["known"]:
        lab = d["known"].split("  ")[0].strip()
        lab = {"Zagier A (7,2,-8)": "Zagier A", "Zagier B (9,3,27)": "Zagier B",
               "Zagier C (10,3,9)": "Zagier C", "Zagier D (11,3,-1)": "Zagier D",
               "Zagier E (12,4,32)": "Zagier E", "Zagier F (17,6,72)": "Zagier F",
               "Apery (17,5,1) zeta(3)/6": "Apéry", "Domb (10,4,64) 7zeta(3)/24": "Domb",
               "T (12,4,16) 7zeta(3)/32": "T", "AZ (9,3,-27) L(3,chi-3)/3": "AZ(9,3,-27)",
               "AZ (11,5,125)": "AZ(11,5,125)", "AZ (7,3,81)": "AZ(7,3,81)",
               "L(f,2) weight-3 cusp row (Domb curve)": "scan-1 $L(f,2)$ row",
               "binom(2n,n)^2-type (2-term)": "2-term degeneration",
               "convolution of binom(2n,n)^2": "self-convolution",
               "Cooper s10-type?": "(scan-1 §4 row)"}.get(lab, lab)
        per = f"**{lab}**" + (f", {per}" if per else "")
    elif not per:
        per = "unidentified" if d.get("digits", 0) >= 20 else "--"
    sl = ",".join(f"$\\sigma_{{{p}}}{{=}}{int(round(v))}$" for p, v in sorted(d["slopes"].items()) if abs(v) > 0.4) or "--"
    return "| %d | %d | $%d,%d,%d$ | $%s$ | %d | $(%d,%d)$ | %.4f | %s | %s | %s | %s | %s | %s | %s | $%s$ | %s |" % (
        i, d["N"], d["M"], d["chi"], d["w"], etastr(d["divs"], d["etaq_t"]), d["tdeg"],
        d["r"]+1, d["D"], d["lam1"],
        ("%.4f" % d["lam2"]) if d["lam2"] is not None else "--",
        ("%.4f" % d["lam_min"]) if d.get("lam_min") is not None else "--",
        ("%g" % d["c"]) if d["c"] is not None else "--", d["k"], sl,
        ("%.3f" % d["score"]) if d["score"] is not None else "--",
        ("**%.3f**" % d["budget"]) if (d["budget"] is not None and d["budget"] > 0.5)
            else (("%.3f" % d["budget"]) if d["budget"] is not None else "--"),
        ",".join(d["a"][:6]), per)

out = [HDR, SEP] + [line(i, d) for i, d in enumerate(rows)]
open(os.path.join(HERE, "table.md"), "w").write("\n".join(out) + "\n")

HDR2 = ("| # | $N$ | $M,\\chi,w$ | $t$ | $\\deg t$ | rec | $\\lambda_1$ | $\\lambda_2$ | "
        "$c$ | $k$ | $\\sigma_p$ | score | budget | $a_n$ | period |")
SEP2 = "|" + "---|"*15
out2, cnt = [HDR2, SEP2], 0
for i, d in enumerate(rows):
    if d["k"] == 0: continue
    parts = line(i, d).split(" | ")
    del parts[8]                                  # drop lambda_min
    parts[-2] = "$" + ",".join(d["a"][:5]) + "$"  # shorter a_n
    out2.append(" | ".join(parts))
    cnt += 1
    if cnt >= 28: break
open(os.path.join(HERE, "top_table.md"), "w").write("\n".join(out2) + "\n")
print(f"table.md: {len(rows)} rows; top_table.md: {cnt} rows; {sum(1 for v in idt.values() if v)} identifications")
