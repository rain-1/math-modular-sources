\\ 06_sweep.gp -- sweep: is F*theta^3*Delta^r in M_{4+12r}(Gamma_0(28)) (trivial or chi_{-4})?
\\ Unknowns: c(-3), c(0) and the DEGEN indices.  r = 3..6.
default(parisize, 6000000000);
NC = 500;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(2500); KNOWN = vector(2500); DEG = List();
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2]);
    if(d>2500, next);
    if(type(v)=="t_STR", listput(DEG,d), CD[d]=v; KNOWN[d]=1)); }
PREC = NC+90;
TH  = Ser(vector(PREC, i, my(n=i-1); if(n==0,1, if(issquare(n), 2, 0))), 'q, PREC);
DE  = 'q*eta('q+O('q^PREC))^24 + O('q^PREC);
T3 = TH^3;
S0 = sum(d=1, NC+30, if(KNOWN[d], CD[d], 0)*'q^d) + O('q^(NC+31));
UN = List(); listput(UN,-3); listput(UN,0);
{ for(i=1,#DEG, if(DEG[i]<=NC+30, listput(UN,DEG[i]))); }
UN = Vec(UN);
print("unknown indices: ", UN);
{ for(r=3,6,
   my(w = 4+12*r, PW = T3*DE^r, RHS = S0*PW, nb, nu=#UN, A, Y, V, M, B);
   for(ci=1,1,
     my(chi = if(ci==1, 1, -4));
     M = if(ci==1, mfinit([28,w],4), mfinit([28,w,-4],4));
     B = mfcoefs(M, NC);
     nb = matsize(B)[2];
     A = matrix(NC+1, nb+nu);
     for(i=1,NC+1, for(j=1,nb, A[i,j]=B[i,j]));
     for(j=1,nu, my(S = 'q^UN[j]*PW + O('q^(NC+1)));
        for(i=1,NC+1, A[i,nb+j] = -polcoeff(S, i-1)));
     Y = vector(NC+1, i, polcoeff(RHS, i-1));
     V = matinverseimage(A, Y~);
     print("r=", r, " weight=", w, " chi=", chi, " dim=", nb, " -> ", if(#V==0, "NO SOLUTION", concat("SOLUTION  c(-3)=", concat(Str(V[nb+1]), concat(" c(0)=", Str(V[nb+2])))))));
  );
}
quit;
