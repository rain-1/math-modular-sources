default(parisizemax, 6000000000);
LOG = "/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_logfit.log";
W(s) = write(LOG, s);
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row1_rec.txt");
Qp = QROW1; RR = 6; QC = vector(RR+1, i, Vecrev(Qp[i]));
default(realprecision, 60);
ev(i,n) = { my(v=QC[i+1], s=0.0, t=1.0); for(j=1,#v, s+=v[j]*t; t*=n); s };
NB = 1000000;
p8 = vector(RR, i, 1.0/8^i);
SS = [2000, 6000, 20000, 60000, 200000, 600000, 1000000];
runf(seed) = { my(cur = vector(NB+1), sm); cur[seed+1] = 1.0/8^seed;
  for(n = seed+1, NB, sm = 0.0;
    for(i = 1, min(RR,n), if(cur[n-i+1]!=0.0, sm += ev(i,n)*p8[i]*cur[n-i+1]));
    cur[n+1] = -sm/(n^4*1.0)); cur};
gettime(); HA = runf(0);
W(Str("=== Row 1 log-fit, N = ", NB, " (A built ", gettime(), " ms) ==="));
W(Str("  (A_n/8^n)*n/(log n)^2 at n=",SS,": ", vector(#SS,i,HA[SS[i]+1]*SS[i]/log(SS[i]*1.0)^2)));
CN = ["B","C","D","E","F"];
for(s=1,5,
  my(h = runf(s), rv, mt, sol);
  rv = vector(#SS, i, h[SS[i]+1]/HA[SS[i]+1]); h = 0;
  W(Str("  ", CN[s], "  r_n = ", rv));
  for(J=2,6,
    my(idx = vector(J+1, i, #SS - J + i));
    mt = matrix(J+1, J+1, a, b, 1.0/log(SS[idx[a]]*1.0)^(b-1));
    sol = matsolve(mt, vector(J+1, a, rv[idx[a]])~);
    W(Str("      fit with ", J+1, " terms in 1/log n -> L = ", sol[1]))));
W("DONE"); quit;
