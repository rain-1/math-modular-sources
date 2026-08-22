\\ 13_baselines.gp -- the two source rows on their OWN hosts (un-Hadamarded),
\\ measured in the same way, as the comparison baseline.
default(parisizemax, 12000000000);
\p 4000
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "rows_raw.gp"));
read(concat(outdir, "rows_full.gp"));
MZ = 1200;
zud(M) = { my(Q = vector(M+1), P = vector(M+1), u, v, w);
  Q[1]=1; Q[2]=7/4; P[1]=0; P[2]=13/8;
  for(m=1,M-1, u=(2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1);
    v=3520*m^6+5632*m^5+2064*m^4-384*m^3-156*m^2+16*m+7;
    w=(2*m-1)^2*(2*m)^2*(20*m^2+32*m+13);
    Q[m+2]=(v*Q[m+1]+w*Q[m])/u; P[m+2]=(v*P[m+1]+w*P[m])/u); [Q,P]; }
ZZ = zud(MZ); ZQ = ZZ[1]; ZP = ZZ[2];
DL = vector(2*MZ+1); DL[1]=1; for(k=1,2*MZ, DL[k+1]=lcm(DL[k],k));

print("=== Zudilin row in its own index m (no 3-section) ===");
print(" m   v2(den Q_m)/m   v2(den P_m)/m   sharp k (D_{2m-1}^k)  Pmax/m   log|Q_m|/m  log|Q_mG-P_m|/m");
{ for(t = 1, 6, my(m = 200*t, q = ZQ[m+1], p = ZP[m+1], dq = denominator(q), dp = denominator(p),
                  od, kk = 0, fa, big = 0);
  od = dp >> valuation(dp,2);
  fa = factor(od); if(matsize(fa)[1] > 0, big = fa[matsize(fa)[1],1]);
  while(kk < 10 && DL[2*m]^kk % od != 0, kk++);
  printf("%4d  %10.4f  %10.4f   %2d   %8.4f   %9.5f  %9.5f\n",
    m, -valuation(dq,2)/m, -valuation(dp,2)/m, kk, 1.0*big/m,
    log(abs(q))/m, log(abs(q*Catalan-p))/m)); }
{ printf(" targets: log|Q|/m -> 5log(phi) = %.5f ; singular points phi^-5 = %.6f, -phi^5 = %.6f\n",
    5*log((1+sqrt(5))/2), ((1+sqrt(5))/2)^(-5), -((1+sqrt(5))/2)^5); }

print("");
print("=== Nesterenko (4,7) row in its own index n ===");
print(" n   v2(den a^N)/n  v2(den b^N)/n  sharp k (D_{6n}^k)  Pmax/n   log|a^N|/n  log|a^N G-b^N|/n");
DT = vector(6*400+1); DT[1]=1; for(k=1,6*400, DT[k+1]=lcm(DT[k],k));
{ for(t = 1, 4, my(n = 100*t, a = AN[n], b = BN[n], da = denominator(a), db = denominator(b),
                  od, kk = 0, fa, big = 0);
  od = db >> valuation(db,2);
  fa = factor(od); if(matsize(fa)[1] > 0, big = fa[matsize(fa)[1],1]);
  while(kk < 10 && DT[6*n+1]^kk % od != 0, kk++);
  printf("%4d  %10.4f  %10.4f   %2d   %8.4f   %9.5f  %9.5f\n",
    n, -valuation(da,2)/n, -valuation(db,2)/n, kk, 1.0*big/n,
    log(abs(a))/n, log(abs(a*Catalan-b))/n)); }

print("");
print("=== Zudilin 3-section (the aligned row used on the Hadamard host) ===");
{ for(t = 1, 4, my(n = 100*t, a = AZ[n], b = BZ[n], da = denominator(a), db = denominator(b),
                  od, kk = 0, fa, big = 0);
  od = db >> valuation(db,2);
  fa = factor(od); if(matsize(fa)[1] > 0, big = fa[matsize(fa)[1],1]);
  while(kk < 10 && DT[6*n+1]^kk % od != 0, kk++);
  printf("%4d  %10.4f  %10.4f   %2d   %8.4f   %9.5f  %9.5f\n",
    n, -valuation(da,2)/n, -valuation(db,2)/n, kk, 1.0*big/n,
    log(abs(a))/n, log(abs(a*Catalan-b))/n)); }
\q
