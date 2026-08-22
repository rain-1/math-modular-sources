\\ 14_pairs.gp -- the adjacent-pair family (j0,j0+1) of POSITIVITY_PROGRAM.md sec 4.3
\\ read as two rows in n, and the CDT data of their Hadamard host.
\\ Fixed rule (no selection bias): m = 3n, j0 = round(0.30 m).
default(parisizemax, 12000000000);
\p 3000
\\ general Catalan-world Beukers moment (copied from lattice/positivity/rows_pos.gp;
\\ that file cannot be read() because its own default(parisizemax,...) aborts the read).
\\ mom(m,j) = I(a,b,a,b,a+1) with a=m+j, b=m; returns [A,B] with mom = A*G - B.
{
mom(m, jj) =
 my(NN = m, TOP = m+jj,
    cc = m!*prod(i=0,m-2,3/2+i)/TOP!,
    num = cc*prod(k=0,TOP-1,('t-k)), tail = vector(NN+2), pre = vector(NN+2),
    MAT, rhs, sol);
 pre[1]=1; for(j=1,NN, pre[j+1]=pre[j]*('t+j));
 tail[NN+2]=1; forstep(j=NN,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 my(bs=vector(2*(NN+1)), lin=vector(NN+2));
 lin[1]=1; for(j=0,NN, lin[j+2]=lin[j+1]*('t+j+1/2));
 for(j=0,NN, bs[2*j+1]=pre[j+1]*lin[j+2]*tail[j+2];
             bs[2*j+2]=pre[j+1]*lin[j+1]*tail[j+2]);
 my(dg=2*NN+1);
 MAT = matrix(dg+1, 2*(NN+1), r, c, polcoeff(bs[c], r-1, 't));
 rhs = vectorv(dg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(MAT, rhs);
 my(A1=vector(NN+1,j,sol[2*(j-1)+1]), A2=vector(NN+1,j,sol[2*(j-1)+2]));
 my(BB = sum(j=0,NN, A2[j+1]),
    CC = sum(j=1,NN, sum(v=0,j-1, (v!/prod(i=0,v-1,3/2+i))*(A1[j+1]+A2[j+1]/(v+1/2)))));
 [4*BB, CC];
}
NTOP = 24;
DT = vector(6*NTOP+2); DT[1] = 1; for(k = 1, 6*NTOP+1, DT[k+1] = lcm(DT[k], k));

print("n  j0  |  log|A1|/n  log|M1|/n | log|A2|/n  log|M2|/n | v2(A1)/n  v2(B1)/n  k1  Pmax1/n | v2(w)/n  log|w|/n  kW");
{ for(n = 2, NTOP,
  my(m = 3*n, j0 = round(0.30*m), r1 = mom(m, j0), r2 = mom(m, j0+1),
     A1 = r1[1], B1 = r1[2], A2 = r2[1], B2 = r2[2],
     M1 = A1*Catalan - B1, M2 = A2*Catalan - B2,
     wn = A1*B2 - A2*B1, dA1 = denominator(A1), dB1 = denominator(B1),
     od, kk = 0, fa, big = 0, dw, odw, kw = 0, vh);
  od = dB1 >> valuation(dB1,2);
  fa = factor(od); if(matsize(fa)[1] > 0, big = fa[matsize(fa)[1],1]);
  while(kk < 12 && DT[6*n+1]^kk % od != 0, kk++);
  dw = denominator(wn); odw = dw >> valuation(dw,2);
  while(kw < 12 && DT[6*n+1]^kw % odw != 0, kw++);
  vh = valuation(numerator(wn),2) - valuation(dw,2);
  printf("%2d %3d | %9.4f %9.4f | %9.4f %9.4f | %8.3f %8.3f  %2d %7.3f | %8.3f %9.4f  kW=%d\n",
    n, j0, log(abs(A1))/n, log(abs(M1))/n, log(abs(A2))/n, log(abs(M2))/n,
    (valuation(numerator(A1),2)-valuation(dA1,2))/n,
    (valuation(numerator(B1),2)-valuation(dB1,2))/n, kk, 1.0*big/n,
    1.0*vh/n, log(abs(wn))/n, kw)); }
\q
