"""Task 3(i): the exact order of the y-ODE satisfied by CDT's G = Sym^+ H,
i.e. how many derivatives of the conditional generator are Q(y)-independent."""
import t3_orbit as T   # reuses G, d_, Int, rank_mod  (prints its own tables first)
p = T.p

def find_relation(fs, dmax=30, L=None):
    """Smallest d such that sum_i P_i(y) f_i = 0 with deg P_i <= d has a nonzero soln."""
    Ly = len(fs[0])
    for d in range(0, dmax+1):
        nun = len(fs)*(d+1)
        L = min(Ly-2, nun+25)
        if L <= nun: return None, None
        rows = []
        for f in fs:
            f = (list(f)+[0]*L)[:L]
            for j in range(d+1): rows.append([0]*j+f[:L-j])
        rk = T.rank_mod(rows, L)
        if rk < nun: return d, (nun-rk)
    return None, None

print()
print("=" * 78)
print("MINIMAL y-ODE ORDER for CDT's conditional generator G  (series to y^%d)" % T.Py)
ONE = [1]+[0]*(T.Py-1)
Gd = [T.G]
for _ in range(9): Gd.append(T.d_(Gd[-1]))
for N in range(1, 9):
    fs = [ONE]+Gd[:N+1]
    d, nul = find_relation(fs)
    print(f"  order N={N}: {{1,G,...,G^({N})}}  ->  " +
          (f"RELATION at deg <= {d} (nullity {nul})" if d is not None else "no relation found in range"))
