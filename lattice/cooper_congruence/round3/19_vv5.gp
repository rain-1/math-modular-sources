\\ 19_vv5.gp -- the vector-valued S-law in the representation forced by the signature.
\\ Disc form (Z/14, Q(x)=x^2/28), signature 7 mod 8 (lattice of signature (1,2)).
\\ For weight 5/2:  2k = 5 = sig + 2 (mod 4)  =>  F_{-beta} = -F_beta  (antisymmetric),
\\ and rho(T) e_beta = e(Q(beta)) e_beta, so the exponents are d/28 with d = beta^2 (28).
\\ Data: c(d) computed with the class bt(d), bt^2 = -3d (28); relabel beta = 3*bt (mod 14),
\\ then beta^2 = 9 bt^2 = -27 d = d (28).
\\ Law tested:  H_b(-1/tau) = C (tau/i)^{5/2} sum_{g=1}^{6} sin(pi b g/7) H_g(tau),  |C| = 1/sqrt(14).
default(realprecision, 60);
DMAX = 2500;
PMAX = 60;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(DMAX); LB = vector(DMAX);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    if(d>DMAX, next);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    if(bt==0 || bt==7, next);
    LB[d] = (3*bt)%14; CD[d] = if(type(v)=="t_STR", 0, v)); }
print("d=1,4,8,9,16,25 -> beta ", vector(6,i,LB[[1,4,8,9,16,25][i]]), "  c ", vector(6,i,CD[[1,4,8,9,16,25][i]]));
{ ampl(beta, lb) = if(lb==0, 0, if(lb==beta, 1, if((14-lb)%14==beta, -1, 0))); }
NEG = List();
{ forstep(d=-1,-PMAX,-1, my(ok=0); for(b=1,6, if((d-b^2)%28==0, ok=b; break)); if(ok, listput(NEG,[d,ok]))); }
NEG = Vec(NEG);
print("admissible negative d (with beta): ", NEG);
{ Hknown(beta,tau) = my(s=0., w=exp(2*Pi*I*tau/28), e);
   for(d=1,DMAX, if(CD[d]!=0, e=ampl(beta,LB[d]); if(e!=0, s += e*CD[d]*w^d))); s; }
{ Hpol(beta,tau,j) = my(w=exp(2*Pi*I*tau/28), e=ampl(beta,NEG[j][2])); if(e==0, 0., e*w^(NEG[j][1])); }
{ SS(b,g) = sin(Pi*b*g/7); }
TAUS = [0.3+0.6*I, -0.2+0.55*I, 0.41+0.72*I, 0.13+0.83*I, -0.37+0.64*I, 0.07+0.91*I, 0.22+0.5*I, -0.44+0.77*I, 0.36+1.05*I, -0.11+0.68*I];
{ run(k) =
  my(nn=#NEG, nu=2*nn+1, nr=6*#TAUS, M=matrix(nr,nu), Y=vector(nr), V, res, C, ok);
  for(ti=1,#TAUS, my(tau=TAUS[ti], t2=-1/tau, pw=exp(k*log(tau/I)));
    for(b=1,6, my(row=6*(ti-1)+b, sa);
      sa = sum(g=1,6, SS(b,g)*Hknown(g,tau));
      for(j=1,nn, M[row,j] = Hpol(b,t2,j));
      M[row,nn+1] = -pw*sa;
      for(j=1,nn, M[row,nn+1+j] = -pw*sum(g=1,6, SS(b,g)*Hpol(g,tau,j)));
      Y[row] = -Hknown(b,t2)));
  if(nr<nu, print("not enough equations: nr=",nr," nu=",nu); return);
  V = matsolve(M[1..nu,], Y[1..nu]~);
  res = M*V - Y~;
  C = V[nn+1];
  ok = vector(nn, j, abs(V[nn+1+j] - C*V[j]));
  print("k=",k,"  maxres=", vecmax(vector(nr,i,abs(res[i]))), "  |C|=", abs(C), "  1/sqrt14=",1/sqrt(14), "  max|w-Cp|=", vecmax(ok));
  if(vecmax(ok) < 1e-15, print("   polar coefficients: ", vector(nn,j,[NEG[j][1], V[j]])));
}
foreach([5/2,3/2,1/2,7/2], k, run(k));
quit;
