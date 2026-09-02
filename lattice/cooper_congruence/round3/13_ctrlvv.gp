\\ 13_ctrlvv.gp -- CONTROL of the vector-valued machinery on a KNOWN object:
\\ the Kohnen cusp form h in S^+_{5/2}(Gamma_0(28)) (symmetric components).
\\ G_b(tau) = sum_{n = b^2 (28)} a(n) q^{n/28};  test G_b(-1/tau) = C (tau/i)^{5/2} sum_g m_g cos(2 pi b g/14) G_g(tau).
default(parisize, 2000000000);
default(realprecision, 60);
NN = 3000;
S = mfinit([28,5/2],1); K = mfkohnenbasis(S); h = mflinear(S, K[,1]);
A = mfcoefs(h, NN);
print("a(1..40): ", vector(40,i,A[i+1]));
{ G(b,tau) = my(s=0., w=exp(2*Pi*I*tau/28)); for(n=1,NN, if((n-b^2)%28==0, s += A[n+1]*w^n)); s; }
MG = [1,2,2,2,2,2,2,1];   \\ multiplicity of gamma=0..7
{ run(tau) = my(t2=-1/tau, pw=exp(5/2*log(tau/I)), L, R, C);
  L = vector(8, i, G(i-1,t2));
  R = vector(8, i, my(b=i-1); pw*sum(j=1,8, my(g=j-1); MG[j]*cos(2*Pi*b*g/14)*G(g,tau)));
  C = L[1]/R[1];
  print("tau=",tau,"  C=",C, "  |C|=",abs(C));
  print("   ratios L/R: ", vector(8,i,if(abs(R[i])<1e-40,"0/0",L[i]/R[i]/C)));
}
run(0.3+0.6*I);
run(-0.2+0.55*I);
quit;
