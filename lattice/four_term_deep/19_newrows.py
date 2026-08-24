#!/usr/bin/env python3
"""19_newrows.py -- analyse a batch of scan output, keep the CANDIDATE rows, and
write their Apery limits in the "label value" format the Catalan battery reads."""
import sys, os, json, importlib.util
HERE=os.path.dirname(os.path.abspath(__file__))
spec=importlib.util.spec_from_file_location("an",os.path.join(HERE,"06_analyse_deep.py"))
an=importlib.util.module_from_spec(spec); spec.loader.exec_module(an)
outjson, outxi = sys.argv[1], sys.argv[2]
paths = sys.argv[3:]
recs = an.main(paths, outjson, full=True)
n=0
with open(outxi,'w') as f:
    for r in recs:
        if r.get('verdict')=='CANDIDATE' and r.get('xi'):
            lab=('x'+'_'.join(map(str,r['cls']))+'__'+'_'.join(map(str,r['row']))).replace('-','m')
            f.write(lab+' '+r['xi']+'\n'); n+=1
print("limits written:", n, file=sys.stderr)
