import json,sys,os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
d=json.load(open('eta_pairs.json'))
os.makedirs('lev',exist_ok=True)
for N,v in d.items():
    if not v['F'] or not v['t']: continue
    with open(f'lev/L{N}.gp','w') as f:
        f.write(f"NLEV={N};\n")
        f.write("DIVS="+str(v['divs']).replace(' ','')+";\n")
        f.write("TS="+str([x['r'] for x in v['t']]).replace(' ','')+";\n")
        f.write("TDEG="+str([x['deg'] for x in v['t']]).replace(' ','')+";\n")
        f.write("FS="+str(v['F']).replace(' ','')+";\n")
    print(N, len(v['t']), len(v['F']))
