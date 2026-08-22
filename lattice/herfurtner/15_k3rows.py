#!/usr/bin/env python3
"""Emit 15_k3run.gp: the K3-window (deg J <= 24) J-map test on

  (a) every row of HERFURTNER_CLASSIFICATION.md sec. 6.1 whose deg <= 12 verdict
      was NEGATIVE  (read off out/jall.log, so the selection is not hand-made);
  (b) the root rows among those (Beukers'/sqrt-Apery, sqrt-T, sqrt-Domb,
      sqrt-AZ(7,3,81), sqrt-s7) -- a subset of (a), flagged in the label;
  (c) three Kodaira-INadmissible controls, expected negative by Theorem H2:
      sqrt-s10, sqrt-s18 (class (8;3,5), delta_inf = 1/4) and the
      NONCONGRUENCE_SCAN.md sec. 4.4 row (class (3;2,5), rho = 7/6);
  (d) the six known POSITIVE rows, as implementation controls: the answer must
      come back with the same (h, deg J) as in out/jall.log.

Run:  python3 15_k3rows.py && gp -q -s 8000000000 15_k3run.gp </dev/null >out/k3jall.log
"""
import os

HERE = os.path.dirname(os.path.abspath(__file__))

NAMES = {
    (1, 0, 0, 7, 2, -8):     "Zagier A",
    (1, 0, 0, 9, 3, 27):     "Zagier B",
    (1, 0, 0, 10, 3, 9):     "Zagier C",
    (1, 0, 0, 11, 3, -1):    "Zagier D",
    (1, 0, 0, 12, 4, 32):    "Zagier E",
    (1, 0, 0, 17, 6, 72):    "Zagier F",
    (2, 1, 1, 20, 2, 16):    "sqrt(Domb)                [root row]",
    (2, 1, 1, 24, 2, 4):     "sqrt(T)                   [root row]",
    (2, 1, 1, 56, 6, 324):   "sqrt(AZ(7,3,81))          [root row]",
    (2, 1, 1, 72, 6, -108):  "sqrt(AZ(9,3,-27))         [root row]",
    (2, 1, 1, 88, 10, 500):  "sqrt(AZ(11,5,125))        [root row]",
    (2, 1, 1, 136, 10, 4):   "Beukers = sqrt(Apery)     [root row]",
    (2, 1, 3, 0, -4, -64):   "A=0 row, class (2;1,3) -",
    (2, 1, 3, 0, 4, -64):    "A=0 row, class (2;1,3) +",
    (3, 1, 1, 117, 21, 441): "NEW row on I1 I7 II II",
    (3, 1, 1, 180, 24, -72): "new row (3;1,1), not a surface",
    (3, 1, 2, 26, 2, -3):    "sqrt(s7) Cooper           [root row]",
    (3, 1, 2, 40, 4, 48):    "new row (3;1,2), not a surface",
    (3, 2, 4, 0, -3, -81):   "A=0 row, class (3;2,4) -",
    (3, 2, 4, 0, 3, -81):    "A=0 row, class (3;2,4) +",
    (4, 1, 3, 28, 2, -8):    "new row (4;1,3), not a surface",
    (4, 1, 3, 48, 4, 32):    "new row (4;1,3), not a surface",
    (4, 1, 3, 68, 6, 72):    "new row (4;1,3), not a surface",
    (4, 1, 3, 72, 6, 108):   "NEW row on I3 III III III",
    (4, 1, 3, 80, 6, 36):    "new row (4;1,3), not a surface",
    (4, 1, 3, 88, 6, -4):    "level-5 Fricke row",
    (4, 3, 5, 0, -4, -256):  "A=0 row, class (4;3,5) -",
    (4, 3, 5, 0, 4, -256):   "A=0 row, class (4;3,5) +",
}

# (c) the Kodaira-inadmissible controls.  P(n) = A(n^2 + (2M-j1-j2)/(2M) n) + B,
#     Q(n) = C (Mn-j1)(Mn-j2).
#  sqrt(s10) : P = 24n^2+12n+2,  Q = -4(8n-3)(8n-5)     SPORADIC_SCAN2 sec. 8
#  sqrt(s18) : P = 56n^2+28n+6,  Q = 12(8n-3)(8n-5)
#  NCS 4.4   : P = 18n^2-3n-6,   Q =  -3(3n-2)(3n-5)    (alpha = 18)
CONTROLS_NEG = [
    ((8, 3, 5, 24, 2, -4),  "sqrt(s10)   [not Kodaira: delta_inf = 1/4]"),
    ((8, 3, 5, 56, 6, 12),  "sqrt(s18)   [not Kodaira: delta_inf = 1/4]"),
    ((3, 2, 5, 18, -6, -3), "NCS 4.4 arctan row, alpha=18  [rho = 7/6]"),
]

# (d) the known positives, as implementation controls
CONTROLS_POS = [
    (1, 0, 0, 9, 3, 27),
    (1, 0, 0, 11, 3, -1),
    (2, 1, 1, 88, 10, 500),
    (2, 1, 1, 72, 6, -108),
    (3, 1, 1, 117, 21, 441),
    (4, 1, 3, 72, 6, 108),
]


def read_jall():
    neg, pos = [], []
    with open(os.path.join(HERE, "out", "jall.log")) as fh:
        for line in fh:
            f = line.split()
            if len(f) != 9:
                continue
            key = tuple(int(z) for z in f[:6])
            (neg if f[6] == "NO" else pos).append(key)
    return neg, pos


def main():
    neg, pos = read_jall()
    rows = []
    for k in CONTROLS_POS:
        assert k in pos, k
        rows.append((k, "CONTROL+ " + NAMES[k]))
    for k in neg:
        rows.append((k, "sec6.1-NEG " + NAMES.get(k, "?")))
    for k, nm in CONTROLS_NEG:
        rows.append((k, "CONTROL- " + nm))

    out = os.path.join(HERE, "15_k3run.gp")
    with open(out, "w") as fh:
        fh.write('read("%s");\n\n' % os.path.join(HERE, "15_k3jtest.gp"))
        fh.write('print("K3 window: dmax=", DMAX, " hmax=", HMAX,'
                 ' " terms=", NTERM, " #gamma=", 4*#GBASE, " p=", KPRIME);\n')
        fh.write('print("label | M j1 j2 A B C | verdict h degJ gamma used extra | ms");\n')
        fh.write('tst(nm, mm, j1, j2, A, B, C) =\n'
                 '{ my(r, ms);\n'
                 '  gettime(); r = k3test(mm, j1, j2, A, B, C); ms = gettime();\n'
                 '  if(r == 0,\n'
                 '     printf("%-42s | %d %d %d %d %d %d | NO  -  -  -  -  - | %d\\n",'
                 ' nm, mm, j1, j2, A, B, C, ms),\n'
                 '     printf("%-42s | %d %d %d %d %d %d | YES %d %d %s %d %d | %d\\n",'
                 ' nm, mm, j1, j2, A, B, C, r[1], r[2], r[3], r[4], r[5], ms));\n'
                 '}\n\n')
        for k, nm in rows:
            fh.write('tst("%s", %d, %d, %d, %d, %d, %d);\n' % ((nm,) + k))
        fh.write("quit;\n")
    print("rows written to 15_k3run.gp:", len(rows))


if __name__ == "__main__":
    main()
