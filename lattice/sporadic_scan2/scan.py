#!/usr/bin/env python3
"""Main scan: for each parameter t and weight w, sweep a box in the integral
lattice of weight-w forms on Gamma_1(N) and keep every F whose t-expansion
a_n satisfies a linear recurrence of order w+1 with polynomial coefficients of
degree <= DMAX (tested by rank deficiency over F_p, batched with numpy).
"""
import sys, os, json, time, argparse
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from qser import divisors, etaquo_series

P = 2147483647           # 2^31-1
NS = 56                  # number of t-expansion coefficients used in the scan
DMAX = 6                 # max degree of the polynomial coefficients

def load_prep(N, w):
    fn = os.path.join(HERE, "prep", f"pr_{N}_{w}.txt")
    if not os.path.exists(fn): return None
    with open(fn) as f:
        lines = [l.strip() for l in f if l.strip()]
    h = lines[0].split()
    dim, rank = int(h[2]), int(h[3])
    if rank < 0: return None
    F0 = np.array(json.loads(lines[1].replace(" ", "")), dtype=object)
    K = [np.array(json.loads(l.replace(" ", "")), dtype=object) for l in lines[2:]]
    return F0, K

def tpow_inv(T, ns, p):
    """Pm[n] = coefficients of t^n ; return Pinv with a_vec = F_vec @ Pinv."""
    Pm = np.zeros((ns, ns), dtype=np.int64)
    cur = np.zeros(ns, dtype=np.int64); cur[0] = 1
    Pm[0] = cur
    for n in range(1, ns):
        new = np.zeros(ns, dtype=np.int64)
        nz = np.nonzero(cur)[0]
        for i in nz:
            new[i:] = (new[i:] + cur[i]*T[:ns-i]) % p
        cur = new
        Pm[n] = cur
    inv = np.zeros((ns, ns), dtype=np.int64)
    for n in range(ns-1, -1, -1):
        inv[n, n] = 1
        for m in range(n+1, ns):
            inv[n, m] = (-(Pm[n, n+1:m+1] @ inv[n+1:m+1, m])) % p
    return Pm, inv

def batch_rank_deficient(M, p, U):
    B, R, C = M.shape
    used = np.zeros((B, R), dtype=bool)
    rank = np.zeros(B, dtype=np.int64)
    ar = np.arange(B)
    rows = np.arange(R)
    for col in range(C):
        colv = np.where(used, 0, M[:, :, col])
        nz = colv != 0
        has = nz.any(axis=1)
        if not has.any():
            break
        pidx = np.argmax(nz, axis=1)
        prow = M[ar, pidx]                      # (B, C)
        piv = prow[:, col]                      # (B,)
        fac = M[:, :, col]                      # (B, R)
        mask = (~used) & (rows[None, :] != pidx[:, None]) & has[:, None]
        upd = (M * piv[:, None, None] - prow[:, None, :]*fac[:, :, None]) % p
        M = np.where(mask[:, :, None], upd, M)
        used[ar, pidx] |= has
        rank += has
        if (rank >= C).all():
            return np.zeros(B, dtype=bool)
    return rank < C

def box_iter(K, B):
    """m-vectors in {-B..B}^K, ordered by increasing sup-norm."""
    import itertools
    vals = [0]
    for b in range(1, B+1): vals += [b, -b]
    for m in itertools.product(vals, repeat=K):
        yield m

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--degmax", type=int, default=2)
    ap.add_argument("--degmin", type=int, default=1)
    ap.add_argument("--kmax", type=int, default=5)
    ap.add_argument("--box", type=int, default=2)
    ap.add_argument("--weights", type=str, default="1,2,3")
    ap.add_argument("--shard", type=int, default=0)
    ap.add_argument("--nshard", type=int, default=1)
    ap.add_argument("--out", type=str, default="hits")
    ap.add_argument("--batch", type=int, default=3000)
    ap.add_argument("--rmax", type=int, default=4)
    args = ap.parse_args()

    tl = json.load(open(os.path.join(HERE, "t_list.json")))
    tl = [x for x in tl if args.degmin <= x["deg"] <= args.degmax]
    weights = [int(x) for x in args.weights.split(",")]
    outfn = os.path.join(HERE, f"{args.out}_{args.shard}.jsonl")
    fout = open(outfn, "a")
    t0 = time.time(); ntested = 0; nhit = 0
    prep_cache = {}
    for ti, tt in enumerate(tl):
        if ti % args.nshard != args.shard: continue
        N = tt["N"]
        T = np.zeros(NS, dtype=np.int64)
        ser = tt["ser"]
        for n in range(1, NS): T[n] = ser[n-1] % P
        Pm, Pinv = tpow_inv(T, NS, P)
        for w in weights:
            key = (N, w)
            if key not in prep_cache:
                prep_cache[key] = load_prep(N, w)
            pr = prep_cache[key]
            if pr is None: continue
            F0, K = pr
            F0m = np.array([int(x) % P for x in F0[:NS]], dtype=np.int64)
            Km = np.array([[int(x) % P for x in k[:NS]] for k in K], dtype=np.int64) if K else np.zeros((0, NS), dtype=np.int64)
            a0 = (F0m @ Pinv) % P
            AK = (Km @ Pinv) % P if len(Km) else np.zeros((0, NS), dtype=np.int64)
            kk = min(len(K), args.kmax)
            r = args.rmax
            DEG = w + 2
            U = (r+1)*(DEG+1)
            NEQ = U + 8
            assert NEQ + r < NS
            npow = np.zeros((NEQ, DEG+1), dtype=np.int64)
            for n in range(1, NEQ+1):
                v = 1
                for j in range(DEG+1):
                    npow[n-1, j] = v; v = v*n % P
            batchM = []; batchm = []
            def flush():
                nonlocal nhit
                if not batchM: return
                arr = np.stack(batchM)                     # (B, NS)
                Bn = arr.shape[0]
                shift = np.stack([arr[:, 1+i:1+i+NEQ] for i in range(r+1)], axis=2)  # (B,NEQ,r+1)
                M = (shift[:, :, :, None] * npow[None, :, None, :]) % P
                M = M.reshape(Bn, NEQ, U)
                defc = batch_rank_deficient(M, P, U)
                for idx in np.nonzero(defc)[0]:
                    fout.write(json.dumps({"ti": ti, "N": N, "r": tt["r"], "D": tt["D"],
                                           "deg": tt["deg"], "w": w,
                                           "m": list(batchm[idx])}) + "\n")
                    nhit += 1
                fout.flush()
                batchM.clear(); batchm.clear()
            for m in box_iter(kk, args.box):
                a = a0.copy()
                for i, mi in enumerate(m):
                    if mi: a = (a + mi*AK[i]) % P
                batchM.append(a); batchm.append(m); ntested += 1
                if len(batchM) >= args.batch: flush()
            flush()
        if ti % 20 == 0:
            print(f"[{time.time()-t0:.0f}s] t {ti}/{len(tl)} N={N} tested={ntested} hits={nhit}", flush=True)
    print(f"DONE tested={ntested} hits={nhit} time={time.time()-t0:.0f}s", flush=True)

if __name__ == "__main__":
    main()
