"""Turn an ident.out line into a readable period expression."""
import re
from fractions import Fraction

NICE = {"zeta(2)": r"\zeta(2)", "zeta(3)": r"\zeta(3)", "zeta(4)": r"\zeta(4)",
        "G=L(2,chi-4)": "G", "L(2,chi-3)": r"L(2,\chi_{-3})", "L(3,chi-3)": r"L(3,\chi_{-3})",
        "L(3,chi-4)": r"L(3,\chi_{-4})", "L(2,chi5)": r"L(2,\chi_5)", "L(3,chi5)": r"L(3,\chi_5)",
        "L(2,chi-7)": r"L(2,\chi_{-7})", "L(2,chi-8)": r"L(2,\chi_{-8})",
        "L(2,chi8)": r"L(2,\chi_8)", "L(2,chi12)": r"L(2,\chi_{12})",
        "L(4,chi-3)": r"L(4,\chi_{-3})", "Pi^3": r"\pi^3", "Pi^2 log2": r"\pi^2\log2",
        "log(2)^3": r"\log^32", "1": "1"}
CUSP = re.compile(r"L\(f_(\d+)w(\d+)c(\d+)n(\d+)e(\d+),(\d+)\)")

def nm(t):
    t = t.strip()
    if t in NICE: return NICE[t]
    m = CUSP.fullmatch(t)
    if m:
        L, k, c, n, e, s = m.groups()
        tag = "" if (n == "1" and e == "1") else f"^{{[{n}.{e}]}}"
        return rf"L(f_{{{L},{k},\chi_{{{c}}}}}{tag},{s})"
    return t

def fmt(msg):
    msg = msg.strip()
    if not msg or msg.startswith("UNIDENT"): return ""
    m = re.fullmatch(r"RATIONAL (-?\d+)/(-?\d+)", msg)
    if m:
        q = Fraction(int(m.group(1)), int(m.group(2)))
        return f"rational ${q.numerator}/{q.denominator}$" if q.denominator != 1 else f"rational ${q.numerator}$"
    if not msg.startswith("x = ("): return msg
    body = msg[len("x = ("):]
    k = body.rfind(")/")
    if k < 0: return msg
    inside, dens = body[:k], body[k+2:]
    try: den = Fraction(int(dens))
    except ValueError: return msg
    parts = []
    chunks = inside.split(" + ")
    try:
        const = Fraction(int(chunks[0]))
    except ValueError:
        return msg
    if const: parts.append((const/den, "1"))
    for ch in chunks[1:]:
        if "*" not in ch: return msg
        co, tn = ch.split("*", 1)
        parts.append((Fraction(int(co))/den, nm(tn)))
    out = []
    for q, tag in parts:
        if q == 0: continue
        if tag == "1":
            out.append(f"{q.numerator}/{q.denominator}" if q.denominator != 1 else f"{q.numerator}")
            continue
        if q == 1: out.append(tag)
        elif q == -1: out.append(f"-{tag}")
        elif q.denominator == 1: out.append(rf"{q.numerator}\,{tag}")
        else:
            sg = "-" if q < 0 else ""
            out.append(rf"{sg}\tfrac{{{abs(q.numerator)}}}{{{q.denominator}}}{tag}")
    s2 = " + ".join(out).replace("+ -", "- ")
    return f"${s2}$" if s2 else ""
