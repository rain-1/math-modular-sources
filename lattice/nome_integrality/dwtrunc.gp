default(parisizemax,8000000000);
N=2600;
A=vector(N+2); A[1]=1; A[2]=5;
{for(n=1,N, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
\\ truncation form:  F_{s+1}(t) F_{s-1}(t^p) = F_s(t) F_s(t^p)  mod p^s   (F_0=1)
Ftr(s,p,K) = sum(n=0,min(p^s-1,K), A[n+1]*t^n) + O(t^(K+1));
sub_tp(f,p,K) = sum(n=0,K\p, polcoeff(f,n)*t^(n*p)) + O(t^(K+1));
{for(pi=1,4, my(p=[2,3,5,7][pi], smax, K);
  smax = 1; while(p^(smax+1) <= N, smax++);
  for(s=1,smax, my(K, L, Rr, D, bad);
    K = min(N, p^(s+1)-1);
    L  = Ftr(s+1,p,K)*sub_tp(Ftr(s-1,p,K),p,K);
    Rr = Ftr(s,p,K)*sub_tp(Ftr(s,p,K),p,K);
    D = L - Rr; bad=0;
    for(i=0,K, if(polcoeff(D,i)%(p^s)!=0, bad++));
    print("Dwork trunc p=",p," s=",s," order<=t^",K,"  violations=",bad)));}
quit;
