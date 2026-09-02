"""Enumerate all intermediate groups Gamma_H = {gamma in Gamma_0(N): d mod N in H},
H a subgroup of (Z/N)^* containing -1, for N <= 60.  Compute mu (index in PSL_2(Z)),
nu_2, nu_3, nu_infty, genus, and the number of special points nu_2+nu_3+nu_infty.

Cosets Gamma_H \ SL_2(Z) <-> primitive (c,d) in (Z/N)^2 modulo (c,d)~(hc,hd), h in H.
Right multiplication acts on the bottom row: (c,d)S=(d,-c), (c,d)T=(c,c+d),
(c,d)ST=(d,d-c) with ST of order 3 in PSL_2.
Output: 01_groups.out  (all genus-0 groups, flagged with #special).
"""
import itertools, json
from math import gcd

def units(N):
    return [a for a in range(N) if gcd(a, N) == 1]

def subgroups_containing_minus1(N):
    """All subgroups H <= (Z/N)^* with -1 in H (as frozensets)."""
    U = units(N)
    if N <= 2:
        return [frozenset(U)]
    # generate all subgroups by closure of subsets of generators (small N: brute force
    # over cyclic subgroups then closure).  |U| <= 20 for N<=60 mostly; do a BFS.
    Uset = set(U)
    subs = set()
    # start from <-1>
    def close(gens):
        H = {1 % N}
        frontier = [1 % N]
        while frontier:
            x = frontier.pop()
            for g in gens:
                y = (x * g) % N
                if y not in H:
                    H.add(y); frontier.append(y)
        return frozenset(H)
    base = close([(-1) % N])
    subs.add(base)
    changed = True
    while changed:
        changed = False
        for H in list(subs):
            for u in U:
                if u in H: continue
                H2 = close(list(H) + [u])
                if H2 not in subs:
                    subs.add(H2); changed = True
    return sorted(subs, key=lambda s: (len(s), sorted(s)))

def cosets(N, H):
    prim = [(c, d) for c in range(N) for d in range(N) if gcd(gcd(c, d), N) == 1]
    rep = {}
    orbits = []
    for p in prim:
        if p in rep: continue
        idx = len(orbits)
        orb = set()
        for h in H:
            orb.add(((h * p[0]) % N, (h * p[1]) % N))
        for q in orb:
            rep[q] = idx
        orbits.append(sorted(orb))
    return rep, orbits

def analyse(N, H):
    rep, orbits = cosets(N, H)
    mu = len(orbits)
    def act(p, M):
        c, d = p
        a, b, cc, dd = M
        return ((c * a + d * cc) % N, (c * b + d * dd) % N)
    S = (0, -1, 1, 0)
    T = (1, 1, 0, 1)
    ST = (0, -1, 1, 1)
    nu2 = 0; nu3 = 0
    for orb in orbits:
        p = orb[0]
        if rep[act(p, S)] == rep[p]: nu2 += 1
        if rep[act(p, ST)] == rep[p]: nu3 += 1
    # T-orbits
    seen = set(); nuinf = 0
    for i, orb in enumerate(orbits):
        if i in seen: continue
        nuinf += 1
        j = i; p = orbits[i][0]
        while True:
            p = act(p, T); j = rep[p]
            if j in seen: break
            seen.add(j)
            if j == i: break
        seen.add(i)
    g2 = 1 + mu / 12 - nu2 / 4 - nu3 / 3 - nuinf / 2
    return mu, nu2, nu3, nuinf, g2

def char_orders(N, H):
    """orders and parities of Dirichlet characters mod N trivial on H:
    i.e. characters of (Z/N)^*/H."""
    U = units(N)
    # quotient group structure: brute force character enumeration via the group (Z/N)^*/H
    cos = {}
    reps = []
    for u in U:
        key = frozenset((u * h) % N for h in H)
        if key not in cos:
            cos[key] = len(reps); reps.append(u)
    n = len(reps)
    return n, reps

if __name__ == '__main__':
    out = []
    for N in range(1, 61):
        for H in subgroups_containing_minus1(N):
            mu, nu2, nu3, nuinf, g = analyse(N, H)
            nsp = nu2 + nu3 + nuinf
            if abs(g) < 1e-9:
                nq, reps = char_orders(N, H)
                out.append(dict(N=N, H=sorted(H), index=nq, mu=mu, nu2=nu2, nu3=nu3,
                                nuinf=nuinf, nspecial=nsp, genus=0))
    with open('01_groups.json', 'w') as f:
        json.dump(out, f, indent=1)
    with open('01_groups.out', 'w') as f:
        f.write("# genus-0 Gamma_H, N<=60.  [(Z/N)^*:H] = index (= #available nebentypus chars)\n")
        f.write("%4s %6s %28s %5s %4s %4s %5s %6s\n" % ("N", "idx", "H", "mu", "nu2", "nu3", "cusps", "#spec"))
        for r in out:
            f.write("%4d %6d %28s %5d %4d %4d %5d %6d%s\n" % (
                r['N'], r['index'], str(r['H'])[:28], r['mu'], r['nu2'], r['nu3'],
                r['nuinf'], r['nspecial'], "   <== FOUR-POINT" if r['nspecial'] == 4 else ""))
    four = [r for r in out if r['nspecial'] == 4]
    print("genus-0 Gamma_H total:", len(out), "  four-point:", len(four))
    for r in four:
        print(r['N'], "H=", r['H'], "idx", r['index'], "mu", r['mu'],
              "nu2,nu3,cusps", r['nu2'], r['nu3'], r['nuinf'])
