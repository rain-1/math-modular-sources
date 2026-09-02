"""Task 3(ii): the ODE order of the X_1(5) Sym^2 row, hence the number of
Q(y)-independent derivatives of its conditional generator after the descent."""
from fractions import Fraction as F
N = 300
a = [F(1), F(3)]
for n in range(1, N):
    a.append((F(11*n*n+11*n+3)*a[n] + F(n*n)*a[n-1])/F((n+1)**2))
assert all(x.denominator == 1 for x in a[:80])
A = [sum(a[i]*a[n-i] for i in range(n+1)) for n in range(N)]      # Sym^2 row: [t^n] F(t)^2
print("Zagier D  a_n :", [int(x) for x in a[:7]])
print("Sym^2     A_n :", [int(x) for x in A[:7]], " (CDT_FINDER: 1,6,47,408,3745,35598,346583)")

p = (1 << 61) - 1
def modp(c): return [((v.numerator % p)*pow(v.denominator % p, p-2, p)) % p for v in c]
def rank_mod(rows, ncol):
    rows = [r[:] for r in rows]; r = 0
    for col in range(ncol):
        piv = next((i for i in range(r, len(rows)) if rows[i][col]), None)
        if piv is None: continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][col], p-2, p); rows[r] = [x*inv % p for x in rows[r]]
        for i in range(len(rows)):
            if i != r and rows[i][col]:
                f0 = rows[i][col]; rows[i] = [(x-f0*y) % p for x, y in zip(rows[i], rows[r])]
        r += 1
        if r == len(rows): break
    return r
def d_(f): return [(n+1)*f[n+1] % p for n in range(len(f)-1)]

for nm, ser in (('weight-1 row F(t)', a), ('Sym^2 row F(t)^2', A)):
    fs = [modp(ser)]
    for _ in range(6): fs.append(d_(fs[-1]))
    print(f"  {nm}: minimal ODE order (homogeneous, poly coefficients):")
    for Nord in range(1, 6):
        found = None
        for d in range(0, 25):
            nun = (Nord+1)*(d+1); L = min(len(fs[Nord])-2, nun+20)
            if L <= nun: break
            rows = []
            for f in fs[:Nord+1]:
                f2 = (list(f)+[0]*L)[:L]
                for j in range(d+1): rows.append([0]*j+f2[:L-j])
            if rank_mod(rows, L) < nun: found = d; break
        print(f"     order {Nord}: " + (f"RELATION at deg <= {found}" if found is not None else "none"))
        if found is not None: break
