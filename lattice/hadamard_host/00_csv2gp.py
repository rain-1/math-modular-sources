#!/usr/bin/env python3
"""Convert the published exact integer rows to PARI/GP vectors."""
import csv, os
base = "/home/ubuntu/code/math-modular-sources/catalan-2-row-denominators/"
out  = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/"
os.makedirs(out, exist_ok=True)
z = list(csv.DictReader(open(base+"zudilin_rows.csv")))
n = list(csv.DictReader(open(base+"nesterenko_rows.csv")))
with open(out+"csvrows.gp","w") as f:
    f.write("XCSV = [" + ",".join(r["X_n"] for r in z) + "];\n")
    f.write("YCSV = [" + ",".join(r["Y_n"] for r in z) + "];\n")
    f.write("VCSV = [" + ",".join(r["V_n"] for r in n) + "];\n")
    f.write("UCSV = [" + ",".join(r["U_n"] for r in n) + "];\n")
print("wrote", out+"csvrows.gp", len(z), len(n))
