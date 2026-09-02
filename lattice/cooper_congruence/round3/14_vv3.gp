\\ 14_vv3.gp -- vector-valued S-law, D-indexing (Jacobi weight 3 index 7 / antisymmetric
\\ weight-5/2 for the Weil rep of Z/14).  Components beta=1..6, exponent D/28 with D = -beta^2 (28);
\\ the data give D = 3d, a(beta,D) = +-c(d).  Unknown polar coefficients at every admissible D<0
\\ down to -PMAX, unknown constant C.  Law: H_b(-1/tau) = C (tau/i)^k sum_g sin(pi b g/7) H_g(tau).
default(realprecision, 70);
DMAX = 2500;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(3*DMAX); LB = vector(3*DMAX);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    if(d>DMAX, next);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    if(bt==0 || bt==7, next);
    LB[3*d] = bt; CD[3*d] = if(type(v)=="t_STR", 0, v)); }
print("D=3,12,24,27,48,75 -> beta ", vector(6,i,LB[[3,12,24,27,48,75][i]]), "  c ", vector(6,i,CD[[3,12,24,27,48,75][i]]));
{ ampl(beta, lb) = if(lb==0, 0, if(lb==beta, 1, if((14-lb)%14==beta, -1, 0))); }
PMAX = 60;
NEG = List();
{ forstep(D=-1,-PMAX,-1, my(ok=0); for(b=1,6, if((D+b^2)%28==0, ok=b; break)); if(ok, listput(NEG,[D,ok]))); }
NEG = Vec(NEG);
print("admissible negative D (with beta): ", NEG);
{ Hknown(beta,tau) = my(s=0., w=exp(2*Pi*I*tau/28), e);
   for(D=1,3*DMAX, if(CD[D]!=0, e=ampl(beta,LB[D]); if(e!=0, s += e*CD[D]*w^D))); s; }
{ Hpol(beta,tau,j) = my(w=exp(2*Pi*I*tau/28), e=ampl(beta,NEG[j][2])); if(e==0, 0., e*w^(NEG[j][1])); }
{ SS(b,g) = sin(Pi*b*g/7); }
TAUS = [0.3+0.6*I, -0.2+0.55*I, 0.41+0.72*I, 0.13+0.83*I, -0.37+0.64*I, 0.07+0.91*I, 0.22+0.5*I];
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
  print("k=",k,"  maxres=", vecmax(vector(nr,i,abs(res[i]))), "  |C|=", abs(C), "  max|w-Cp|=", vecmax(ok));
  if(vecmax(ok) < 1e-20, print("   polar coefficients p: ", vector(nn,j,[NEG[j][1], V[j]])));
}
foreach([5/2,3/2,7/2,1/2], k, run(k));
quit;
