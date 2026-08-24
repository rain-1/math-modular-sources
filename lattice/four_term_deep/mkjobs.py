#!/usr/bin/env python3
"""Emit the job queue for a deep four-term scan, ordered by a-block (so that a
prefix of the queue is a complete scan of a low-a box: 'prioritise by a then d')."""
import sys
BLOCKS = [1,30,70,120,180,250,330,420,520,640,800,1000,1300,1700,2200]
def blocks(amax):
    out=[]
    for i in range(len(BLOCKS)-1):
        lo, hi = BLOCKS[i]+ (0 if i==0 else 1), BLOCKS[i+1]
        if lo>amax: break
        out.append((lo, min(hi,amax)))
    return out

mode   = sys.argv[1]          # DG | DW | DWB | MIX
cfile  = sys.argv[2]
amax   = int(sys.argv[3])
CMAX,DMAX,FMAX,GMAX = (int(x) for x in sys.argv[4:8])
N      = int(sys.argv[8]) if len(sys.argv)>8 else 24
cls=[l.split() for l in open(cfile) if l.strip()]
jobs=[]
for (lo,hi) in blocks(amax):
    for c in cls:
        tag = "_".join(c).replace("-","m")
        if mode=='MIX':
            # c = RPN RPD RRN RRD M J1 J2 ; shard by r-range instead of a
            jobs.append(f"./03_fmix {' '.join(c)} {lo} {hi} {amax} {CMAX} {DMAX} {FMAX} {GMAX} {N} > out/{mode}_{tag}_rp{lo}-{hi}.txt 2>logs/{mode}_{tag}_rp{lo}-{hi}.log")
            jobs.append(f"./03_fmix {' '.join(c)} {-hi} {-lo} {amax} {CMAX} {DMAX} {FMAX} {GMAX} {N} > out/{mode}_{tag}_rn{lo}-{hi}.txt 2>logs/{mode}_{tag}_rn{lo}-{hi}.log")
        else:
            D = DMAX if mode=='DG' else -1
            G = GMAX if mode=='DG' else -1
            jobs.append(f"./02_fscan {' '.join(c)} {lo} {hi} {CMAX} {D} {FMAX} {G} {N} > out/{mode}_{tag}_a{lo}-{hi}.txt 2>logs/{mode}_{tag}_a{lo}-{hi}.log")
print("\n".join(jobs))
