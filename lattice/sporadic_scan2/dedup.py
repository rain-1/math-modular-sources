#!/usr/bin/env python3
"""Deduplicate verified rows, flag the known census rows, rank by budget/score,
and emit limits.txt for the PARI identification pass."""
import sys, os, json, math, argparse

HERE = os.path.dirname(os.path.abspath(__file__))

KNOWN = {
 (1,2,10,56,346):      "Zagier A (7,2,-8)  zeta(2)/4",
 (1,3,9,21,9):         "Zagier B (9,3,27)",
 (1,3,15,93,639):      "Zagier C (10,3,9)  L(2,chi-3)/2",
 (1,3,19,147,1251):    "Zagier D (11,3,-1) zeta(2)/5",
 (1,4,20,112,676):     "Zagier E (12,4,32) G/2",
 (1,6,42,312,2394):    "Zagier F (17,6,72)",
 (1,5,73,1445,33001):  "Apery (17,5,1) zeta(3)/6",
 (1,4,28,256,2716):    "Domb (10,4,64) 7zeta(3)/24",
 (1,4,40,544,8536):    "T (12,4,16) 7zeta(3)/32",
 (1,3,27,309,4059):    "AZ (9,3,-27) L(3,chi-3)/3",
 (1,5,35,275,2275):    "AZ (11,5,125)",
 (1,3,9,3,-279):       "AZ (7,3,81)",
 (1,2,12,104,1078):    "L(f,2) weight-3 cusp row (Domb curve)",
 (1,4,4,16,36):        "binom(2n,n)^2-type (2-term)",
 (1,8,88,1088,14296):  "convolution of binom(2n,n)^2",
 (1,2,10,60,386):      "Cooper s10-type?",
}

def canon(a):
    a = [int(x) for x in a]
    b = [((-1)**n)*a[n] for n in range(len(a))]
    for x, y in zip(a, b):
        if x != y:
            return tuple(a) if x > y else tuple(b)
    return tuple(a)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rows", nargs="+", required=True)
    ap.add_argument("--out", default="table.json")
    args = ap.parse_args()
    seen = {}
    for fn in args.rows:
        for line in open(fn):
            d = json.loads(line)
            key = canon(d["a"][:8])
            cur = seen.get(key)
            if cur is None or (d["N"], d["M"], d["w"]) < (cur["N"], cur["M"], cur["w"]):
                if cur is not None:
                    d["also"] = cur.get("also", 0) + 1
                seen[key] = d
            else:
                cur["also"] = cur.get("also", 0) + 1
    rows = list(seen.values())
    for d in rows:
        d["known"] = KNOWN.get(canon(d["a"][:5]), "")
        d["canon"] = list(canon(d["a"][:8]))
    rows.sort(key=lambda d: -(d["budget"] if d["budget"] is not None else -99))
    json.dump(rows, open(os.path.join(HERE, args.out), "w"), indent=0)
    with open(os.path.join(HERE, "limits.txt"), "w") as f:
        for i, d in enumerate(rows):
            if d.get("limit") and d.get("digits", 0) >= 20 and not d.get("cplx"):
                f.write(f"{i} {d['M']} {d['w']} {min(d['digits'],200)} {d['limit']}\n")
    print(f"{len(rows)} distinct rows; {sum(1 for d in rows if d['known'])} known-flagged")

if __name__ == "__main__":
    main()
