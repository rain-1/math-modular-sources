default(realprecision, 3000);
G = Catalan;
NBQ(i,K,Q) = my(c=0); for(h=1,K, if((2*i+2*h+1)%Q==0, c++)); c;
vpfact(n,p) = my(v=0,m=n); while(m>0, m=m\p; v+=m); v;
tails(Gv,n) = my(T=vector(n),part=0); for(m=1,n, part += (-1)^(m-1)/(2*(m-1)+1)^2; T[m] = (-1)^m*(Gv-part)); T;
mQ(Q,B,S,N,A) = my(cA,cost,cls,marg,v); cA = sum(k=1,S, (A[k]+2*B)\Q); cost = vector(N, ii, 2*NBQ(ii-1,B,Q) - NBQ(ii-1,S,Q) - 2*(Q<=2*(ii-1)+1) - ((ii-1)\Q) - ((N-ii)\Q)); cls = vector(Q, r, []); for(ii=1,N, cls[((ii-1)%Q)+1] = concat(cls[((ii-1)%Q)+1], [cost[ii]])); marg = []; for(r=1,Q, v=vecsort(cls[r]); for(k=1,#v, marg = concat(marg, [v[k]+2*(k-1)]))); marg = vecsort(marg); cA + sum(k=1,S, marg[k]);
test(B) =
{
  my(S,N,D,Tm,um,PP,R,best,bestA,A,detR,logF,logPi,logH,M,X,Ap,mm,Q,delta0);
  S = B\20; N = 2*B+S+3; D = 2*B;
  Tm = tails(G, N+S+1);
  um = vector(N+S+1, m, Tm[m]/(2*m+1));
  PP = vector(N, i, prod(h=1,B,(2*(h+i-1)+1)^2));
  R = matrix(S+3,S, a, j, sum(i=0,a-1+2*B, (-1)^i*binomial(a-1+2*B,i)*PP[i+1]*um[i+j]));
  best = -1; bestA = 0;
  forsubset([S+3,S], sub, my(d=abs(matdet(matrix(S,S,k,j,R[sub[k],j])))); if(d>best, best=d; bestA=sub));
  A = vector(S, k, bestA[k]-1);
  detR = best;
  logF = sum(r=0,D-1, lngamma(r+1));
  logPi = sum(i=1,N, log(PP[i]));
  logH = 0; M = 0;
  forprime(p=3, 6*B+20,
    Ap = sum(i=1,N, valuation(PP[i],p)) - sum(r=0,D-1, vpfact(r,p));
    mm = 0; Q = p;
    while(Q <= 6*B+20, mm += mQ(Q,B,S,N,A); Q *= p);
    M += mm*log(p);
    logH += max(Ap - mm, 0)*log(p));
  X = log(detR) + logF - logPi + logH;
  delta0 = 0.00966242652523235;
  printf("B=%d S=%d N=%d  log|detR|=%.3f  M_odd=%.3f  v2(F)log2=%.3f  X_B=%.3f  -delta0*B^2=%.3f  X/B^2=%.5f\n", B, S, N, log(detR), M, sum(r=0,D-1,vpfact(r,2))*log(2), X, -delta0*B^2, X/B^2);
}
for(k=1,6, test(20*k));
test(160); test(200);
quit;
