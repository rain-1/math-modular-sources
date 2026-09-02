"""For every four-point genus-zero curve X_H (from 01_groups.py) list every nebentypus
character that can occur on a group Gamma_H' with the SAME image in PSL_2(Z)
(i.e. H'*{+-1} = H), together with its order and parity.  A source Phi with
non-real Fourier coefficients requires a nebentypus of order >= 3.

Output 03_chars.out
"""
from math import gcd
import cmath, json

def units(N): return [a for a in range(N) if gcd(a,N)==1]

def subgroups(N):
    U = units(N); out=set()
    def close(gens):
        H={1%N}; fr=[1%N]
        while fr:
            x=fr.pop()
            for g in gens:
                y=(x*g)%N
                if y not in H: H.add(y); fr.append(y)
        return frozenset(H)
    out.add(close([]))
    changed=True
    while changed:
        changed=False
        for H in list(out):
            for u in U:
                if u in H: continue
                H2=close(list(H)+[u])
                if H2 not in out: out.add(H2); changed=True
    return out

def chars_trivial_on(N,H):
    """characters of (Z/N)^* trivial on H, as functions given by values on units,
    returned as (order, parity, valuestring)."""
    U = units(N)
    # build the abelian group (Z/N)^*/H explicitly and its dual by brute force over
    # all maps: enumerate characters via cyclic decomposition using primitive roots
    # -- simpler: enumerate all functions chi:U->roots of unity that are homomorphisms
    # by taking the quotient group and listing its characters through a Smith-like
    # generator search.
    # quotient elements
    cos={}; reps=[]
    for u in U:
        key=frozenset((u*h)%N for h in H)
        if key not in cos: cos[key]=len(reps); reps.append(u)
    n=len(reps)
    def cls(u): return cos[frozenset((u*h)%N for h in H)]
    # multiplication table -> find generators greedily
    gens=[]; gen_orders=[]; covered={cls(1)}
    for u in reps:
        c=cls(u)
        if c in covered: continue
        # order of u in quotient
        o=1; v=u%N
        while cls(v)!=cls(1): v=(v*u)%N; o+=1
        newcov=set()
        for a in covered:
            ra=reps[a]
            for j in range(o):
                newcov.add(cls((ra*pow(u,j,N))%N))
        if len(newcov)>len(covered):
            gens.append(u); gen_orders.append(o); covered=newcov
        if len(covered)==n: break
    # characters: assign gens[i] -> zeta_{gen_orders[i]}^{k_i}; may over-count if the
    # decomposition is not direct, so filter for well-definedness.
    import itertools
    out=[]
    for ks in itertools.product(*[range(o) for o in gen_orders]):
        # define chi on the subgroup generated; check well-defined by testing on all U
        val={}
        ok=True
        # enumerate all products of generators
        for exps in itertools.product(*[range(o) for o in gen_orders]):
            u=1
            for g,e in zip(gens,exps): u=(u*pow(g,e,N))%N
            c=cls(u)
            z=sum(ks[i]*exps[i]/gen_orders[i] for i in range(len(gens)))%1.0
            if c in val and abs(val[c]-z)>1e-9 and abs(abs(val[c]-z)-1)>1e-9:
                ok=False; break
            val[c]=z
        if not ok or len(val)!=n: continue
        # order & parity
        from fractions import Fraction
        fr=[Fraction(round(v*1_000_000),1_000_000).limit_denominator(1000) for v in val.values()]
        order=1
        for f in fr: order=order*f.denominator//gcd(order,f.denominator)
        m1=(-1)%N
        p=val[cls(m1)]
        parity = 1 if abs(p)<1e-9 or abs(p-1)<1e-9 else (-1 if abs(p-0.5)<1e-9 else 0)
        out.append((order,parity))
    return sorted(set(out))

GR = json.load(open('01_groups.json'))
FOUR = [(g['N'], tuple(g['H']), 'N=%d H=%s  (#special=%d%s)' % (g['N'], g['H'], g['nspecial'],
        ', FOUR-POINT' if g['nspecial']==4 else '')) for g in GR]

lines=[]
lines.append("# nebentypus characters available on EVERY genus-zero Gamma_H curve, N <= 60")
lines.append("# (all H' with H'*{+-1} = Hbar, i.e. same curve; order 3,4,6 = non-real coefficients)")
for N,Hbar,name in FOUR:
    Hbar=frozenset(x%N for x in Hbar)
    lines.append("")
    lines.append(f"{name}   N={N}  Hbar={sorted(Hbar)}")
    allch=set()
    for H in subgroups(N):
        Hb=frozenset(list(H)+[(-x)%N for x in H])
        if Hb!=Hbar: continue
        cs=chars_trivial_on(N,H)
        lines.append(f"   H'={str(sorted(H)):<20} -> characters (order,parity): {cs}")
        allch|=set(cs)
    orders=sorted(set(o for o,_ in allch))
    lines.append(f"   ==> character orders available: {orders};  max order = {max(orders)}"
                 f"{'   <== NON-REAL (order>=3) POSSIBLE' if max(orders)>=3 else '   (all real: order<=2)'}")
open('03_chars.out','w').write("\n".join(lines)+"\n")
print("\n".join(lines))
