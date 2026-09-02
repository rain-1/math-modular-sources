/* 08_al.gp -- Atkin-Lehner action on u for the three rows, numerically (60 digits).
   u_{s7}  = (eta(7t)/eta(t))^4
   u_{s10} = (eta(5t)eta(10t)/(eta(t)eta(2t)))^2
   u_{s18} = (eta(2t)eta(3t)^2 eta(18t)/(eta(t)eta(6t)^2 eta(9t)))^6
   W_Q = [Q a, b; N c, Q d], det Q, Q||N.  u is weight 0.                       */
default(realprecision, 60);
{ RW = [["s7",7,13,49],["s10",10,6,25],["s18",18,14,1]]; }
NAM=["s7","s10","s18"];
{ uu(k,t) =
  if(k==1, return( (eta(7*t,1)/eta(t,1))^4 ));
  if(k==2, return( (eta(5*t,1)*eta(10*t,1)/(eta(t,1)*eta(2*t,1)))^2 ));
  (eta(2*t,1)*eta(3*t,1)^2*eta(18*t,1)/(eta(t,1)*eta(6*t,1)^2*eta(9*t,1)))^6;
}
{ ALmat(N,Q) = my(R=N/Q);
  for(a=1,60, for(d=-60,60, for(b=-60,60, for(c=-60,60,
    if(Q*a*Q*d - b*N*c == Q, return([Q*a,b;N*c,Q*d]))))));
  0; }
print("=== sanity: u(tau) at tau = 0.13+0.37i vs the q-series ===");
{ for(k=1,3, print("  ",NAM[k],": u = ",uu(k,0.13+0.37*I))); }
print();
{ for(k=1,3, my(N=RW[k][2], C=RW[k][4], t=0.13+0.37*I, u0);
  u0 = uu(k,t);
  print("=== ",NAM[k],"  N=",N,"  C=",C);
  fordiv(N, Q, if(Q==1 || gcd(Q,N/Q)!=1, next);
    my(W=ALmat(N,Q), wt, uw);
    wt = (W[1,1]*t+W[1,2])/(W[2,1]*t+W[2,2]);
    uw = uu(k,wt);
    print("   Q=",Q,"  W=",W);
    print("      u|W_Q / u      = ",uw/u0);
    print("      u|W_Q * (C u)  = ",uw*C*u0);
    if(abs(uw/u0-1)<1e-40, print("      ==> u|W_",Q," = u ;  F|W = +F, rho|W = +rho, Phi|W_",Q," = +Phi"));
    if(abs(uw*C*u0-1)<1e-40, print("      ==> u|W_",Q," = 1/(Cu) ; F|W = -F, rho|W = -rho, Phi|W_",Q," = -Phi"));
  )); }
quit;
