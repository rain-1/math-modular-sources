\\ Task 1: the level-60 purified weight-8 Eisenstein source for zeta(7).
\\ Phi = sum_{d|60} c_d E_8(d tau);  L(Phi,s) = P(s) zeta(s) zeta(s-7),
\\ P(s) = sum_d c_d d^{-s}.  Purity <=> P(0)=P(2)=P(4)=P(6)=P(8)=0.
\\ L(Phi,7) = (-P(7)/2) zeta(7)   since zeta(0) = -1/2.
default(realprecision,50);
dl = divisors(60);
P(c,s) = sum(i=1,#dl, c[i]*dl[i]^(-s));
pure(c) = vector(5,j,P(c,2*j-2))==[0,0,0,0,0];
print("divisors of 60: ", dl);

\\ (a) p=5 completion Phi60^(5) = (1 - 5^4 V_5) Phi12
c12 = [1,-572,11583,-36608,46332,-20736]; d12=[1,2,3,4,6,12];
{c5 = vector(#dl, i, my(d=dl[i], v=0);
      for(j=1,6, if(d12[j]==d, v+=c12[j]); if(5*d12[j]==d, v-=5^4*c12[j]));
      v);}
print("Phi60^(5) = ", c5);
print("  = archive vector? ", c5==[1,-572,11583,-36608,-625,46332,357500,-20736,-7239375,22880000,-28957500,12960000]);
print("  purified (kills s=0,2,4,6,8)? ", pure(c5));
print("  -P(7)/2 = ", -P(c5,7)/2, "  vs 6479/54000 = ", 6479/54000);
print("  level-12 parent -P12(7)/2 = ", -sum(j=1,6,c12[j]*d12[j]^(-7))/2, " vs 209/1728=",209/1728);

\\ (b) Atkin-Lehner: E_k(d tau)|W_Q = (Q/a^2)^{k/2} E_k(dQ/a^2 tau), a = Q-part of d
{AL(c,Q) = my(k=8, out=vector(#dl));
  for(i=1,#dl, my(d=dl[i], a=gcd(d,Q^10), dn, f, j=0);
    dn = d*Q/a^2; f = (Q/a^2)^(k/2);
    for(t=1,#dl, if(dl[t]==dn, j=t));
    out[j] = out[j] + f*c[i]);
  out;}
{for(qi=1,7, my(Q=[3,4,5,12,15,20,60][qi], w=AL(c5,Q));
  print("  W_",Q," eigen: +1? ", w==c5, "   -1? ", w==-c5));}

\\ (c) purified space, and its W_60 = +1 part
M = matrix(5,#dl, r,i, dl[i]^(-(2*(r-1))));
{A60 = matrix(#dl,#dl); for(i=1,#dl, my(e=vector(#dl)); e[i]=1;
   my(w=AL(e,60)); for(j=1,#dl, A60[j,i]=w[j]));}
print("dim purified (12-dim oldform space): ", #matker(M));
print("dim purified & W60=+1 : ", #matker(matconcat([M; A60-matid(#dl)])));
print("dim purified & W60=-1 : ", #matker(matconcat([M; A60+matid(#dl)])));

\\ (d) archive basis
B0=[0,0,2673,-73216,398125,-497664,-3840000,13208832,-14478750,5280000,0,0];
B1=[0,33,0,-39424,276250,-375921,-2900625,9165312,-7796250,0,1670625,0];
B2=[1,0,0,-402688,3062500,-4313088,-33280000,101606400,-79633125,0,0,12960000];
{for(i=1,3, my(B=[B0,B1,B2][i]);
  print("B",i-1," pure? ",pure(B), "  W60=+1? ",AL(B,60)==B, "  -P(7)/2 = ", -P(B,7)/2));}
print("L(B0/480,7)/zeta(7) = ", -P(B0,7)/2/480, "  vs 2623/216000 = ", 2623/216000);
print("L(B0/216,7)/zeta(7) = ", -P(B0,7)/2/216);
MB = Mat([B0~,B1~,B2~]);
print("rank(B0,B1,B2) = ", matrank(MB));
sol = matsolve(MB~*MB, MB~*c5~);
print("Phi60^(5) = ", sol~, " . (B0,B1,B2) ? ", MB*sol==c5~);

\\ (e) support-depth filtration
{for(mm=1,7, my(rows=[]);
  for(i=1,#dl, if(dl[i]<mm, rows=concat(rows,[i])));
  my(Mx = if(#rows==0, matconcat([M; A60-matid(#dl)]),
              matconcat([M; A60-matid(#dl); matrix(#rows,#dl,r,i, if(i==rows[r],1,0))])));
  print("dim(purified & W60=+1 & c_d=0 for d<",mm,") = ", #matker(Mx)));}
\q
