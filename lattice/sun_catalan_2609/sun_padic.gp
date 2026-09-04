default(realprecision, 200);
NBQ(i,K,Q) = my(c=0); for(h=1,K, if((2*i+2*h+1)%Q==0, c++)); c;
tails(Gv,n) = my(T=vector(n),part=0); for(m=1,n, part += (-1)^(m-1)/(2*(m-1)+1)^2; T[m] = (-1)^m*(Gv-part)); T;
mQ(Q,B,S,N,A) = my(cA,cost,cls,marg,v); cA = sum(k=1,S, (A[k]+2*B)\Q); cost = vector(N, ii, 2*NBQ(ii-1,B,Q) - NBQ(ii-1,S,Q) - 2*(Q<=2*(ii-1)+1) - ((ii-1)\Q) - ((N-ii)\Q)); cls = vector(Q, r, []); for(ii=1,N, cls[((ii-1)%Q)+1] = concat(cls[((ii-1)%Q)+1], [cost[ii]])); marg = []; for(r=1,Q, v=vecsort(cls[r]); for(k=1,#v, marg = concat(marg, [v[k]+2*(k-1)]))); marg = vecsort(marg); cA + sum(k=1,S, marg[k]);
padtest(B, Gp) =
{
  my(S,N,q,Tm,um,PP,R,A,detR,mm,Q,ok,bad);
  S = B\20; N = 2*B+S+3; q = denominator(Gp);
  Tm = tails(Gp, N+S+1);
  um = vector(N+S+1, m, Tm[m]/(2*m+1));
  PP = vector(N, i, prod(h=1,B,(2*(h+i-1)+1)^2));
  R = matrix(S+3,S, a, j, sum(i=0,a-1+2*B, (-1)^i*binomial(a-1+2*B,i)*PP[i+1]*um[i+j]));
  A = vector(S,k,k+2);
  detR = matdet(matrix(S,S,k,j,R[A[k]+1,j]));
  if(detR==0, print("det zero for this A"); return);
  ok = 1; bad = [];
  forprime(p=3, 6*B+20,
    mm = 0; Q = p;
    while(Q <= 6*B+20, mm += mQ(Q,B,S,N,A); Q *= p);
    if(valuation(q^S*detR, p) < mm, ok = 0; bad = concat(bad, [[p, valuation(q^S*detR,p), mm]])));
  printf("B=%d S=%d log10(q)~%d : lemma holds at all odd p? %d ; violations (p, v_p, claimed bound): %s\n", B, S, round(log(q)/log(10)), ok, bad);
}
padtest(40, bestappr(Catalan, 10^12));
padtest(40, bestappr(Catalan, 10^30));
padtest(60, bestappr(Catalan, 10^20));
padtest(80, 7/8);
quit;
