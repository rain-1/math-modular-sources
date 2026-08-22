import json,os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
d=json.load(open('eta_pairs.json'))
CH=20000   # pairs per job
jobs=[]
for N,v in sorted(d.items(), key=lambda x:int(x[0])):
    nt,nf=len(v['t']),len(v['F'])
    if nt==0 or nf==0: continue
    step=max(1, CH//max(nf,1))
    for lo in range(1,nt+1,step):
        jobs.append((int(N),lo,min(lo+step-1,nt)))
os.makedirs('out/mod',exist_ok=True)
with open('jobs.txt','w') as f:
    for i,(N,lo,hi) in enumerate(jobs):
        f.write(f"{i} {N} {lo} {hi}\n")
print(len(jobs),"jobs")
