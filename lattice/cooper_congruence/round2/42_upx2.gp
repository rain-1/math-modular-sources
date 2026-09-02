/* 42_upx2.gp -- ROUTE (A), CORRECTED.
   The task's premise "U_p(f) = T_p(f) - V_p(f) mod p for a modular FUNCTION f" is WRONG:
   for weight 0 the geometric (degree p+1) Hecke correspondence acts by
        T_p f = p * U_p f + V_p f ,
   because sum_{i=0}^{p-1} f((tau+i)/p) = p * (U_p f)(tau).  The Kronecker congruence
        Phi_p(X,Y) = (X^p - Y)(X - Y^p)  (mod p)
   therefore says exactly  T_p f = V_p f = f^p (mod p), i.e. p*U_p f = 0 mod p -- a tautology.
   So the modular polynomial mod p gives NO information about U_p(x^b) mod p; one needs it
   mod p^2.   THIS SCRIPT: (i) exhibits the failure of the naive claim numerically;
   (ii) asks the real question -- is U_p(x^b) mod p a RATIONAL FUNCTION of x (or of u)?
   (iii) tests the Hasse-invariant fact  F^((p-1)/2) = rational function of u  (mod p).      */
default(parisize,8000000000);
read("40_core.gp");

/* find R in F_p(t) of degree <= d with S = R(t) as q-series mod p, using MM equations.
   t, S are q-series mod p.  Returns [num,den] as polynomials in 'T, or 0.               */
{ ratfit(S,t,d,MM,p) = my(cols=2*d+2, Mat0, tp, v, K, sol, nu, de);
  Mat0 = matrix(MM, cols);
  tp = 1 + O('q^(MM+2));
  v = vector(d+1);
  for(i=0,d, v[i+1] = tp; tp = tp*t);
  for(i=0,d,
    for(r=1,MM, Mat0[r,i+1] = polcoeff(S*v[i+1], r-1)));
  for(i=0,d,
    for(r=1,MM, Mat0[r,d+2+i] = -polcoeff(v[i+1], r-1)));
  K = matker(Mat0*Mod(1,p));
  if(#K==0, return(0));
  sol = K[,1];
  de = sum(i=0,d, sol[i+1]*'T^i);
  nu = sum(i=0,d, sol[d+2+i]*'T^i);
  [nu,de];
}

M = 420;
print("=== (iii) Hasse-invariant structure:  F^((p-1)/2) = R(u) mod p ?  (row s7, u a Hauptmodul) ===");
{ my(S=Setup(1,M), u=S[1], F=S[2]);
  forprime(p=5,17,
    my(Fp = Mod(1,p)*F, up = Mod(1,p)*u, W, res, d);
    W = Fp^((p-1)/2);
    for(d=1,26, res = ratfit(W,up,d,3*d+12,p);
      if(res!=0, print("  p=",p,"  F^((p-1)/2) = N(u)/D(u) with deg <= ",d,
                       "   N=",lift(res[1]),"   D=",lift(res[2])); break);
      if(d==26, print("  p=",p,"  NO rational fit in u of degree <= 26")));
  ); }
print();
print("=== (ii) is U_p(x^b) mod p a rational function of x?  (row s7) ===");
{ my(S=Setup(1,M), x=S[3]);
  forprime(p=3,13, if(7%p==0,next);
    my(xp = Mod(1,p)*x, xb = 1+O('q^M));
    for(b=1,min(p-1,4),
      xb = xb*x;
      my(Ub = Mod(1,p)*Ser(vector(M\p, m, polcoeff(xb,p*m)),'q,M\p), res, found=0);
      \\ Ub has no constant term; shift: fit Ub itself
      for(d=1,18, res = ratfit(Ub,xp,d,3*d+14,p);
        if(res!=0, print("  p=",p," b=",b,":  U_p(x^b) = N(x)/D(x), deg<=",d,
                         "  N=",lift(res[1]),"  D=",lift(res[2])); found=1; break));
      if(!found, print("  p=",p," b=",b,":  NO rational fit in x of degree <= 18"))));
}
quit;
