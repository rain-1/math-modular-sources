\\ 59_w2s10.gp -- the Atkin-Lehner W_2 action on the s10 Heegner set, and the factor 2.
\\ Phi|W_2 = +Phi, f|W_2 = +f, so fhat is W_2-INVARIANT: the natural home of the trace is
\\ X_0(10)/W_2.  Question: does W_2 act freely on {disc -4m^2, 10|A, B=6m (20)}/Gamma_0(10),
\\ and does it preserve chi_{-4} ?  If yes for m odd and not for m even, the factor 2 is explained.
read("50_lib.gp"); read("heeg.gp"); read("maass2.gp");
default(realprecision, 50);
N=10; initfser(2,700);
bet = read("20_beta_s10.txt");
cont(Q) = gcd(gcd(Q[1],Q[2]),Q[3]);
\\ smallest-A Gamma_0(N) Heegner representative of Q
{ heegmin(Q,N,beta,PB) = my(a=Q[1],b=Q[2],c=Q[3],A,B,C,g,s,q,best=0);
  for(r=0,PB, for(p=-PB,PB,
    if(gcd(p,r)!=1, next);
    A = a*p^2+b*p*r+c*r^2;
    if(A<=0 || A%N!=0, next);
    if(best!=0 && A>=best[1], next);
    g = bezout(p,-r); if(g[3]!=1, next);
    s = g[1]; q = g[2];
    B = 2*a*p*q + b*(p*s+q*r) + 2*c*r*s;
    if((B-beta)%(2*N)!=0, next);
    C = (B^2-(b^2-4*a*c))/(4*A);
    best = [A,B,C]));
  best;
}
\\ W_2 = [[6,1],[10,2]], det 2, normalises Gamma_0(10).  Image form of a Heegner form.
{ w2form(Qh, d) = my(A=Qh[1],B=Qh[2], al, al2, A2, B2, C2);
  al = (-B + I*sqrt(-d))/(2*A);
  al2 = (6*al+1)/(10*al+2);
  A2 = round(sqrt(-d)/(2*imag(al2)));
  B2 = round(-2*A2*real(al2));
  C2 = (B2^2-d)/(4*A2);
  if(denominator(C2)!=1, return(0));
  [A2,B2,C2];
}
{ redof(Q) = my(R=qfbred(Qfb(Q[1],Q[2],Q[3]))); my(v=Vec(R)); [v[1],v[2],v[3]]; }
{
for(m=1,20,
  if(m%5==0, next);
  my(d=-4*m^2, bt=(6*m)%20, RF=redforms(d), REP=vector(#RF), CH=vector(#RF), FH=vector(#RF), OM=vector(#RF), img=vector(#RF), T=0., nfix=0, T2=0., seen=vector(#RF));
  for(i=1,#RF,
    REP[i] = heegmin(RF[i],N,bt,60);
    if(REP[i]==0, next);
    CH[i] = genchar(REP[i],-4); OM[i] = omeg(REP[i]);
    FH[i] = fhatR(2, (-REP[i][2]+I*2*m)/(2*REP[i][1]));
    T += CH[i]*FH[i]/OM[i]);
  \\ W_2 action on the index set (via SL_2(Z)-reduced classes)
  for(i=1,#RF,
    if(REP[i]==0, img[i]=0; next);
    my(W=w2form(REP[i],d), r0, j0=0);
    if(W==0, img[i]=0; next);
    r0 = redof(W);
    for(j=1,#RF, if(RF[j]==r0, j0=j; break));
    img[i]=j0);
  \\ transversal sum
  for(i=1,#RF,
    if(REP[i]==0 || seen[i], next);
    seen[i]=1; if(img[i]>0, seen[img[i]]=1);
    if(img[i]==i, nfix++);
    T2 += CH[i]*FH[i]/OM[i]);
  print("m=",m,"  #cls=",#RF,"  W_2-fixed classes=",nfix,
        "   chi/fhat preserved by W_2? ",
        my(ok=1); for(i=1,#RF, if(REP[i]==0||img[i]==0, next);
           if(CH[i]!=CH[img[i]] || abs(FH[i]-FH[img[i]])>1e-20*(1+abs(FH[i])), ok=0)); ok);
  print("     T      = ",T,"    beta/(i T)  = ",bet[m]/(I*T));
  print("     T_W2   = ",T2,"    beta/(i T_W2) = ",bet[m]/(I*T2),"     beta=",bet[m]);
  print("     W_2 permutation: ",img);
);
}
quit;
