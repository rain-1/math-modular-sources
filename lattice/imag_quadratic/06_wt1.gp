\\ 06_wt1.gp -- weight-one rows (the CDT architecture) on the four-point hosts that
\\ admit an odd nebentypus, with every admissible pole placement (hostscan sec.3),
\\ including the placements at a complex-conjugate special point (K imaginary quadratic).
\\ Calibration: Gamma_0(6) -> Zagier A/C/F ; Gamma_0(8) -> Zagier E.
default(seriesprecision,90);
M=90; q='q;
w3 = Mod('t,'t^2+3);

minrec(A,name,lam) = {
  print("== ",name);
  print("   A_0..A_6 = ", vector(7,i,A[i]));
  my(found=0);
  for(T=2,6, for(deg=1,6,
    my(cols=T*(deg+1), rows=40, m);
    if(rows+T+2>#A, next);
    m=matrix(rows,cols);
    for(e=1,rows, my(n=e+1);
      for(j=0,T-1, for(dd=0,deg, m[e,j*(deg+1)+dd+1]=n^dd*A[n+j+1])));
    my(k=matker(m));
    if(#k>0 && !found,
      print("   MINIMAL recurrence: terms=",T," deg=",deg," (kernel dim ",#k,")");
      if(T==3, my(v=k[,1],c0,c1,c2);
        c0=v[deg+1]; c1=v[2*(deg+1)]; c2=v[3*(deg+1)];
        print("   char. poly  ",lift(c2),"*y^2 + ",lift(c1),"*y + ",lift(c0),
              "   roots(y=1/lambda)=",polroots(lift(c2)*'y^2+lift(c1)*'y+lift(c0)));
        print("   expected 1/lambda = ",lam));
      found=1)));
  if(!found, print("   no recurrence with terms<=6, deg<=6"));
}

\\ ---------- calibration: Gamma_0(6), weight 1, chi_-3 ----------
G6=znstar(6,1); mf6=mfinit([6,1,[G6,[1]]],3);
print("dim M_1(Gamma_1(6),chi_-3) = ",mfdim(mf6));
F6 = Ser(mfcoefs(mfbasis(mf6)[1],M),q);
h6 = eta(q^2+O(q^M))*eta(q^6+O(q^M))^5/(eta(q+O(q^M))^5*eta(q^3+O(q^M)))*q;
print("h6 = ",h6+O(q^6),"   (mu: cusp0=0, 1/2=-9, 1/3=-8)");
minrec(Vec(subst(F6,q,serreverse(h6))),                 "Gamma_0(6) pole cusp 0  = Zagier F (lambda=-9,-8)", [-1/9,-1/8]);
minrec(Vec(subst(F6,q,serreverse(h6/(1+9*h6)))),        "Gamma_0(6) pole 1/2     = Zagier C = CDT (lambda=9,1)", [1/9,1]);
minrec(Vec(subst(F6,q,serreverse(h6/(1+8*h6)))),        "Gamma_0(6) pole 1/3     = Zagier A (lambda=8,-1)", [1/8,-1]);

\\ ---------- the candidate: Gamma_0(7)-curve, weight 1, chi_-7 ----------
print();
G7=znstar(7,1); mf7=mfinit([7,1,[G7,[3]]],3);
print("dim M_1(Gamma_1(7),chi_-7) = ",mfdim(mf7));
F7 = Ser(mfcoefs(mfbasis(mf7)[1],M),q);
print("F7 = ",F7+O(q^8));
h7 = q*eta(q^7+O(q^M))^4/eta(q+O(q^M))^4;
minrec(Vec(subst(F7,q,serreverse(h7))), "X_0(7) pole cusp 0  (lambda = -13/2 +- 3sqrt3/2 i, |.|=7 both)", 0);
minrec(Vec(subst(F7,q,serreverse(h7/(1-((-13+3*w3)/2)*h7)))),
       "X_0(7) pole ell_3  K=Q(sqrt-3)  (lambda = (13-3sqrt-3)/2, -3sqrt-3)", [2/(13-3*sqrt(-3)), -1/(3*sqrt(-3))]);

\\ ---------- Gamma_0(9)-curve, weight 1, chi_-3 ----------
print();
G9=znstar(9,1); mf9=mfinit([9,1,[G9,[3]]],3);
print("dim M_1(Gamma_1(9),chi_-3 mod 9) = ",mfdim(mf9));
F9 = Ser(mfcoefs(mfbasis(mf9)[1],M),q);
h9 = q*eta(q^9+O(q^M))^3/eta(q+O(q^M))^3;
minrec(Vec(subst(F9,q,serreverse(h9))), "X_0(9) pole cusp 0 = Zagier B (lambda=-9/2 +- 3sqrt3/2 i)", 0);
minrec(Vec(subst(F9,q,serreverse(h9/(1-((-9-3*w3)/2)*h9)))),
       "X_0(9) pole cusp 1/3  K=Q(sqrt-3) (lambda = (9+3sqrt-3)/2, 3sqrt-3 : EQUAL moduli)", 0);
quit;
