\\ 08_sweep2.gp -- wider sweep.  F = sum_d c(d) q^d with unknowns at d<=0 admissible and at
\\ the DEGEN indices.  Test F*theta^{8-j}*Delta^r in M_{4+12r}(Gamma_0(L)) for half-weight j/2,
\\ j=1,3,5,7; levels L=28,56,112; principal parts PMAX = 3 or 31.
default(parisize, 8000000000);
NC = 400; PREC = NC+120;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(2500); KNOWN = vector(2500); DEG = List();
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2]);
    if(type(v)=="t_STR", listput(DEG,d), CD[d]=v; KNOWN[d]=1)); }
TH = Ser(vector(PREC, i, my(n=i-1); if(n==0,1, if(issquare(n),2,0))), 'q, PREC);
DE = 'q*eta('q+O('q^PREC))^24 + O('q^PREC);
S0 = sum(d=1, NC+60, if(KNOWN[d], CD[d], 0)*'q^d) + O('q^(NC+61));
sq28 = Set(vector(28,i,(i-1)^2 % 28));
{ negadm(PM) = my(L=List()); forstep(d=0,-PM,-1, if(setsearch(sq28, (d%28)), listput(L,d))); Vec(L); }
{ trial(j, r, L, PM) =
   my(w = 4+12*r, PW = TH^(8-j)*DE^r, UN, nu, M, B, nb, A, Y, V, RHS);
   UN = concat(negadm(PM), select(x->x<=NC+60, Vec(DEG)));
   nu = #UN;
   RHS = S0*PW + O('q^(NC+1));
   M = mfinit([L,w],4); B = mfcoefs(M,NC); nb = matsize(B)[2];
   A = matrix(NC+1, nb+nu);
   for(i=1,NC+1, for(jj=1,nb, A[i,jj]=B[i,jj]));
   for(jj=1,nu, my(S = 'q^UN[jj]*PW + O('q^(NC+1)));
      for(i=1,NC+1, A[i,nb+jj] = -polcoeff(S,i-1)));
   Y = vector(NC+1, i, polcoeff(RHS,i-1));
   V = matinverseimage(A, Y~);
   print("j=",j," wt=",j,"/2  r=",r," L=",L," PM=",PM," dim=",nb," nunk=",nu," -> ",
         if(#V==0,"NO", concat("YES  c(-3)=", Str(V[nb+#negadm(PM)]))));
}
{ forstep(j=1,7,2, trial(j,3,28,3)); }
print("--- bigger principal part (PM=31), needs Delta^3*Delta(28tau): use r=3 and L=28 with PM=3 only; instead raise r ---");
{ forstep(j=1,7,2, trial(j,4,28,3)); }
{ forstep(j=1,7,2, trial(j,3,56,3)); }
{ forstep(j=1,7,2, trial(j,3,112,3)); }
quit;
