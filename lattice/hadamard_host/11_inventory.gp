\\ 11_inventory.gp -- the full candidate inventory on the Hadamard host, with
\\ archimedean growth, 2-adic slope, sharp LCM exponent, and singular set.
default(parisizemax, 12000000000);
\p 6000
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "rows_raw.gp"));
read(concat(outdir, "rows_full.gp"));
read(concat(outdir, "had.gp"));
N = #HZZ;
AB = bestappr(Catalan, 10^8); aa = numerator(AB); bb = denominator(AB);
{ printf("test rational a/b = %d/%d   |G-a/b| = %.3e\n", aa, bb, 1.0*abs(Catalan-AB)); }

CZq = vector(N, n, bb*BZ[n] - aa*AZ[n]);     \\ conditional Zudilin form (rational)
CNq = vector(N, n, bb*BN[n] - aa*AN[n]);     \\ conditional Nesterenko form (rational)
CONDZ = vector(N, n, AZ[n]*CNq[n]);          \\ A_Z (*) (bB_N - aA_N)
CONDN = vector(N, n, AN[n]*CZq[n]);          \\ A_N (*) (bB_Z - aA_Z)
DBLq  = vector(N, n, CZq[n]*CNq[n]);         \\ doubly conditional

DT = vector(6*N+1); DT[1] = 1; for(k = 1, 6*N, DT[k+1] = lcm(DT[k], k));
Dn(n) = DT[6*n+1];

NMS = ["1          ","A_Z(*)A_N  ","A_Z(*)B_N  ","A_N(*)B_Z  ","B_Z(*)B_N  ","W          ","COND_Z     ","COND_N     ","DBL        "];
SQ  = [vector(N,n,if(n==0,1,0)), HZZ, HZB, HNB, HBB, WW, CONDZ, CONDN, DBLq];

print("");
print("function      growth(n=400)  2-adic slope  sharp k (D_{6n}^k)  b (Pmax/n)");
{ for(s = 2, 9,
    my(c = SQ[s][400], dd = denominator(c), v2 = valuation(dd,2), od, fa, big = 0, kk = 0);
    od = dd >> v2;
    fa = factor(od); if(matsize(fa)[1] > 0, big = fa[matsize(fa)[1],1]);
    kk = 0; while(kk < 10 && Dn(400)^kk % od != 0, kk++);
    printf("%s  %10.5f   %10.5f     %2d              %7.4f\n",
      NMS[s], log(abs(c))/400,
      (valuation(numerator(c),2)-v2)/400, kk, 1.0*big/400)); }

print("");
print("2-adic slope stability (v_2(c_n)/n at n = 100,200,300,400):");
{ for(s = 2, 9,
    printf("%s", NMS[s]);
    for(k = 1, 4, my(n = 100*k, c = SQ[s][n]);
      printf("  %9.4f", (valuation(numerator(c),2)-valuation(denominator(c),2))/n));
    print("")); }

print("");
print("archimedean growth log|c_n|/n at n = 100,200,300,400:");
{ for(s = 2, 9,
    printf("%s", NMS[s]);
    for(k = 1, 4, my(n = 100*k); printf("  %9.5f", log(abs(SQ[s][n]))/n));
    print("")); }
\q
