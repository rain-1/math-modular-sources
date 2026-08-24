#!/usr/bin/env python3
"""06_sizes.py -- the conformal (cusp) size at x=0 of the POST-HYPOTHESIS orbifold
of each mixed-exponent four-term host, by the numerical Fuchsian uniformisation of
05_uniformise.py.

Post-hypothesis geometry: the fold is deleted (H = B - xi A is holomorphic there,
verified in 03_foldreg.py); what remains is P^1 with a cusp at 0, a cusp at
infinity, and the two other finite singular points with their own types (cusp if
the exponent difference is 0, order-2 cone point if it is 1/2).

Reported per row:
  size            the computed cusp size |varphi'(0)|
  ceil16          16 * min |t| over the post-hypothesis set  (rigorous
                  monotonicity ceiling: fill in one of the two points)
  loss            log(ceil16/size), the cost of the fourth special point
  size_pre_ceil   16 * |t_fold|, the ceiling WITHOUT the hypothesis
"""
import os, sys, json, math, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location("uni", os.path.join(HERE, "05_uniformise.py"))
uni = importlib.util.module_from_spec(spec); spec.loader.exec_module(uni)

SPANS = [1.0, 4.0, 16.0, 64.0]


def size_of(pts, ngrid=600, verbose=False):
    best = None
    sc = 1.0 / max(1e-12, abs(complex(pts[0][0])) * abs(complex(pts[1][0])))
    for f in SPANS:
        try:
            r = uni.solve_size(pts, ngrid=ngrid, span=max(f * sc, f * 0.25))
        except Exception as e:
            if verbose:
                print("    span %.4g: %s" % (f * sc, e))
            continue
        if verbose:
            print("    span %.4g: size %.9f (acc %+.8f)" % (max(f*sc, f*0.25), r['size'], r['acc']))
        if best is None or r['size'] > best['size']:
            best = r
    if best is None:
        raise RuntimeError("no size found")
    return best


def main():
    geo = json.load(open(os.path.join(HERE, 'out', 'geometry.json')))
    out = []
    for g in geo:
        pts = g['pts']
        fold = pts[0]
        post = pts[1:]
        cfg = [(complex(p['tre'], p['tim']), 0.0 if p['tag'] == 'cusp' else 0.5) for p in post]
        print("%s  post = %s" % (g['label'],
              ", ".join("%s%+.6fi [%s]" % (("%.6f" % p['tre']), p['tim'], p['tag']) for p in post)))
        t = time.time()
        try:
            r = size_of(cfg, verbose=True)
        except Exception as e:
            print("   FAILED: %s" % e); continue
        loss = math.log(r['ceil16'] / r['size'])
        rec = dict(label=g['label'], size=r['size'], acc=r['acc'], imc=r['imc'],
                   ceil16=r['ceil16'], loss=loss,
                   logsize=math.log(r['size']),
                   pre_ceil=16 * fold['absl'], t_fold=fold['absl'])
        out.append(rec)
        print("   size = %.9f   log = %+.6f   ceiling %.6f (log %+.6f)   loss %.4f nats   [%.0fs]"
              % (r['size'], math.log(r['size']), r['ceil16'], math.log(r['ceil16']), loss, time.time()-t))
    json.dump(out, open(os.path.join(HERE, 'out', 'sizes.json'), 'w'), indent=1)
    print("\nwrote out/sizes.json")


if __name__ == '__main__':
    main()
