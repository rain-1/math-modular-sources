\\ 22_basis.gp -- the integral basis of J^!_{3,7} = M^{!,-}_{5/2}(rho_L), L the Gamma_0(7)
\\ Heegner lattice (disc form Z/14, x^2/28).  For each admissible D0 < 0 we look for the form
\\ g_{D0} with principal part exactly (e_{beta} - e_{-beta}) q^{D0/28} and no other D<0 term,
\\ inside Delta^{-n} J^weak_{3+12n,7}; when it exists we print its coefficients and test
\\ integrality.  Also records which single-slot principal parts are realizable (obstructions).
default(parisize, 12000000000);
read("15_jaclib.gp");
P21 = mkphi21(); P01 = mkphi01(); P12 = mkphi12();
{ Mbasis(w) = my(L=List());
  if(w<0 || w%2, return([]));
  if(w==0, return([one]));
  if(w==2, return([]));
  for(j=0, w\12, my(rem = w-12*j);
    for(b=0,1, my(r2 = rem-6*b);
      if(r2>=0 && r2%4==0, listput(L, DELS^j * E4S^(r2/4) * E6S^b))));
  Vec(L); }
{ Jbasis(n) = my(L=List(), base, MB);
  for(b=0,5,
    base = jmul(P12, jmul(jpow(P01,5-b), jpow(P21,b)));
    MB = Mbasis(4+2*b+12*n);
    for(i=1,#MB, listput(L, jscal(MB[i], base))));
  Vec(L); }
{ dinv(n) = my(D = DELS^n, u = vector(NQ+1), s);
  for(i=0,NQ, u[i+1] = polcoeff(D, 2*(i+n), 'y));
  s = Ser(u, 'y, NQ+1);
  Vec(1/s + O('y^(NQ+1))); }
\\ admissible (D, beta) with beta in 1..7 (beta and 14-beta identified), D = -beta^2 (28)
{ slots(lo, hi) = my(L=List());
  for(D=lo,hi, for(b=1,6, if((D+b^2)%28==0, listput(L,[D,b]); break)));
  Vec(L); }
{ Fc(BS, IU, n, i, D, r) = my(N=(D+r^2)/28);
  if(N+n < 0, return(0));
  sum(j=0, N+n, IU[j+1]*jcoef(BS[i], N+n-j, r)); }

NLEVEL = 1;
BS = Jbasis(NLEVEL); nb = #BS; IU = dinv(NLEVEL);
NEGS = slots(-49-28*NLEVEL, -1);
print("n=", NLEVEL, "  dim=", nb, "  negative slots: ", NEGS);
A = matrix(#NEGS, nb, i, j, Fc(BS, IU, NLEVEL, j, NEGS[i][1], NEGS[i][2]));
print("rank of the principal-part map = ", matrank(A), " (dim = ", nb, ")");
POS = slots(1, 28*(NQ-NLEVEL-2));
print("number of positive slots computed: ", #POS);
{ for(t=1,#NEGS,
   my(e = vector(#NEGS, i, if(i==t,1,0)), V, cf, den);
   V = matinverseimage(A, e~);
   if(#V==0, print("D0=", NEGS[t][1], " (beta=", NEGS[t][2], "): NOT realizable at n=", NLEVEL); next);
   cf = vector(min(24,#POS), i, sum(j=1,nb, V[j]*Fc(BS,IU,NLEVEL,j,POS[i][1],POS[i][2])));
   den = 1; for(i=1,#cf, den = lcm(den, denominator(cf[i])));
   print("D0=", NEGS[t][1], " (beta=", NEGS[t][2], "):  denominator of the first 24 coefficients = ", den);
   print("     c(D) at D=", vector(min(8,#POS),i,POS[i][1]), " : ", vector(min(8,#cf),i,cf[i])));
}
quit;
