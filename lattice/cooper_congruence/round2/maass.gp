\\ maass.gp -- for each Cooper row: the weight -2 weakly holomorphic form f = 1/(xF)
\\ on Gamma_0(N) (q-expansion q^{-1}+O(1), f|W_N = -f, holomorphic on H) and the
\\ Gamma_0(N)-invariant real-analytic weight-0 function
\\        fhat = D f + f/(2 pi y) = -R_{-2} f/(4 pi).
read("lib.gp"); read("e2.gp");
{ E4hol(z,PR) = my(q=exp(2*Pi*I*z), s=1.); for(n=1,PR, s += 240*sigma(n,3)*q^n); s; }
{ E4v(t) = my(R=sl2red(t), z, c, d, PR);
  z=R[1]; c=R[4]; d=R[5];
  PR = ceil(1.2*default(realprecision)/(-log(abs(exp(2*Pi*I*z)))/log(10))) + 40;
  E4hol(z,PR)/(c*t+d)^4;
}
{ E2v(t) = E2star(t) + 3/(Pi*imag(t)); }
\\ eta-quotient u for each row, evaluated at a point
{ uv(k,t) = if(k==1, (eta(7*t,1)/eta(t,1))^4,
             if(k==2, (eta(5*t,1)*eta(10*t,1)/(eta(t,1)*eta(2*t,1)))^2,
                      (eta(2*t,1)*eta(3*t,1)^2*eta(18*t,1)/(eta(t,1)*eta(6*t,1)^2*eta(9*t,1)))^6));
}
\\ F = sum_d lam_d E_2(d tau) ; lam lists (cooper_sources 3.2)
{LAMS = [[[1,-1/6],[7,7/6]],
         [[1,-1/12],[2,-2/12],[5,5/12],[10,10/12]],
         [[1,-1/4],[2,2/4],[3,6/4],[6,-12/4],[9,-9/4],[18,18/4]]];}
{ Fv(k,t) = my(L=LAMS[k], s=0.); for(i=1,#L, s += L[i][2]*(E2star(L[i][1]*t) + 3/(Pi*imag(L[i][1]*t)))); s; }
{ Fmod(k,t) = my(L=LAMS[k], s=0.); for(i=1,#L, s += L[i][2]*E2star(L[i][1]*t)); s; }
{ DFv(k,t) = my(L=LAMS[k], s=0., d); for(i=1,#L, d=L[i][1]; s += L[i][2]*d*(E2v(d*t)^2-E4v(d*t))/12); s; }
\\ closed-form fhat (fails only where F=0)
{ fhatC(k,t) = my(R=ROWS[k], B=R[3], C=R[4], u, g, x, F, Dx, DF, f, Df);
  u = uv(k,t); g = 1+B*u+C*u^2; x = u/g;
  F = Fmod(k,t);
  Dx = F*u*(1-C*u^2)/g^2;
  DF = DFv(k,t);
  f = 1/(x*F);
  Df = -f^2*(F*Dx + x*DF);
  [f, Df + f/(2*Pi*imag(t))];
}
\\ q-series fallback (used at the zeros of F, i.e. at the polar CM points)
{ mkfser(k,NQ) = my(S=Setup(k,NQ), f); f = 1/(S[3]*S[2]); [f, 'q*deriv(f,'q)]; }
{ evser(f,q,TT) = my(v=valuation(f,'q), s=0.); for(n=v,TT, s += polcoeff(f,n)*q^n); s; }
\\ maximise Im over <T, [1,0;kN,1], W_N>, tracking the sign from f|W_N = -f
{ redN(t,N) = my(s=1, n, best, bt, bs, cand, ii);
  for(it=1,200,
    n = round(real(t)); t = t-n;
    best = imag(t); bt = t; bs = s; ii = 0;
    for(kk=-4,4, if(kk==0, next);
      cand = t/(kk*N*t+1);
      if(imag(cand) > best*(1+1e-30), best=imag(cand); bt=cand; bs=s; ii=1));
    cand = -1/(N*t);
    if(imag(cand) > best*(1+1e-30), best=imag(cand); bt=cand; bs=-s; ii=1);
    if(ii==0, break);
    t = bt; s = bs);
  [t,s];
}
{ fhatQ(t,N,FF,TT) = my(R=redN(t,N), z, q, y);
  z=R[1]; y=imag(z); q=exp(2*Pi*I*z);
  R[2]*( evser(FF[2],q,TT) + evser(FF[1],q,TT)/(2*Pi*y) );
}
