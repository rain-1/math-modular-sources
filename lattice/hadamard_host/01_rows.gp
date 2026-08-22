\\ 01_rows.gp -- exact rational rows for the Hadamard host of (Zudilin, Nesterenko)
\\ Normalisation:  az_n := Q_{3n}, bz_n := P_{3n}   (bz/az -> G)
\\                 an_n := 4*B_n,  bn_n := C_n      (bn/an -> G)
default(parisizemax, 12000000000);
\p 4000
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "csvrows.gp"));

NMAX = 400;

zud(M) =
{ my(Q = vector(M+1), P = vector(M+1), aa, bb, cc);
  Q[1] = 1; Q[2] = 7/4; P[1] = 0; P[2] = 13/8;
  for(m = 1, M-1,
    aa = (2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1);
    bb = 3520*m^6 + 5632*m^5 + 2064*m^4 - 384*m^3 - 156*m^2 + 16*m + 7;
    cc = (2*m-1)^2*(2*m)^2*(20*m^2+32*m+13);
    Q[m+2] = (bb*Q[m+1] + cc*Q[m])/aa;
    P[m+2] = (bb*P[m+1] + cc*P[m])/aa);
  [Q, P];
}

nestB(n) =
{ my(s = 0, t);
  for(j = 0, 3*n,
    t = FT[8*n+2*j+1] * FT[j+1] * FT[6*n+1];
    t /= (FT[4*n+1] * FT[4*n+j+1] * FT[3*n-j+1]^2 * FT[2*j+1]^2);
    s += t * 2^(-14*n + 2*j + 1));
  s;
}

print("factorial table ...");
FT = vector(14*NMAX + 3); FT[1] = 1;
for(k = 1, 14*NMAX+2, FT[k+1] = FT[k]*k);

print("Zudilin recurrence to m = ", 3*NMAX, " ...");
ZZ = zud(3*NMAX); ZQ = ZZ[1]; ZP = ZZ[2];
AZ = vector(NMAX, n, ZQ[3*n+1]);
BZ = vector(NMAX, n, ZP[3*n+1]);

print("Nesterenko B_n ...");
BB = vector(NMAX, n, nestB(n));
AN = vector(NMAX, n, 4*BB[n]);

DT = vector(6*NMAX+1); DT[1] = 1;
for(k = 1, 6*NMAX, DT[k+1] = lcm(DT[k], k));
Dn(n) = DT[6*n+1];
ez(n) = my(m = 3*n); min(6*m, 4*m + 3 + logint(2*m-1, 2));

bad = 0;
{ for(n = 1, #XCSV,
  if(XCSV[n] != 2^ez(n)*Dn(n)^2*AZ[n], bad++; print("  X mismatch n=", n));
  if(YCSV[n] != 2^ez(n)*Dn(n)^2*BZ[n], bad++; print("  Y mismatch n=", n));
  if(VCSV[n] != 4^(7*n+1)*Dn(n)^2*BB[n], bad++; print("  V mismatch n=", n))); }
{ print("cross-check vs published integer rows (n<=", #XCSV, "): ",
      if(bad, concat(Str(bad), " MISMATCHES"), "ALL AGREE")); }

BN0 = vector(#UCSV, n, UCSV[n] / (4^(7*n) * Dn(n)^2));

write(concat(outdir, "rows_raw.gp"), "AZ = ", AZ, ";");
write(concat(outdir, "rows_raw.gp"), "BZ = ", BZ, ";");
write(concat(outdir, "rows_raw.gp"), "AN = ", AN, ";");
write(concat(outdir, "rows_raw.gp"), "BN0 = ", BN0, ";");
print("wrote out/rows_raw.gp");

{ for(k = 1, 4, my(n = 100*k);
    printf("n=%3d  log|az|/n=%.6f  log|an|/n=%.6f  log|az G-bz|/n=%.6f\n",
      n, log(abs(AZ[n]))/n, log(abs(AN[n]))/n, log(abs(AZ[n]*Catalan-BZ[n]))/n)); }
{ for(k = 1, 4, my(n = 24*k);
    printf("n=%3d  log|an G-bn|/n=%.6f\n", n, log(abs(AN[n]*Catalan-BN0[n]))/n)); }
{ printf("15log(phi) = %.6f   2logT=%.6f  2logT^- =%.6f\n",
  15*log((1+sqrt(5))/2), 2*log((3303+437*sqrt(57))/144), 2*log((3303-437*sqrt(57))/144)); }
\q
