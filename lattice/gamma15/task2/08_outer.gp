/* 08_outer.gp -- ODE right-hand sides for ALL FOUR real Eisenstein directions
   R1,R2 (outer) and R3,R4 (inner), to locate the pole of the inhomogeneity. */
default(parisizemax, 20000000000);
NA = 90; NQ = NA+3;
leg5(n) = kronecker(n,5);
re4(n) = my(r=n%5); if(r==1, 1, if(r==4, -1, 0));
im4(n) = my(r=n%5); if(r==2, 1, if(r==3, -1, 0));
Ap(n) = sum(k=0, n, binomial(n,k)^2*binomial(n+k,k));
xser(nq) = q*prod(n=1, nq-1, (1 - q^n + O(q^nq))^(5*leg5(n)));
Dinv2(s) = my(m=serprec(s,q)); sum(k=1, m-1, polcoeff(s,k)/k^2*q^k) + O(q^m);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
xs = xser(NQ); Fs = sum(n=0, NQ-1, Ap(n)*xs^n) + O(q^NQ); Av = peel2(Fs, xs, NA, NQ);
mk(f) = sum(n=1, NQ-1, f(n)*q^n) + O(q^NQ);
R1 = 2*mk(n->sumdiv(n,d,re4(d)*d^2));
R2 = -2*mk(n->sumdiv(n,d,im4(d)*d^2));
R3 = 2*mk(n->sumdiv(n,d,re4(n/d)*d^2));
R4 = -2*mk(n->sumdiv(n,d,im4(n/d)*d^2));
Lc(bb, n) = my(b0=bb[n+1]); my(b1=if(n>=1, bb[n], 0)); my(b2=if(n>=2, bb[n-1], 0)); n^2*b0 - (11*n^2-11*n+3)*b1 - (n-1)^2*b2;
rep(Ph, nm) = my(B=peel2(Fs*Dinv2(Ph), xs, NA, NQ)); my(S=sum(i=1,NA, Lc(B,i)*x^(i-1)) + O(x^NA)); my(N=S*(1-11*x-x^2)); print(nm, ":  B_1..B_4 = ", vector(4,i,B[i+1]), "    RHS*(1-11x-x^2) = ", vector(6,i,polcoeff(N,i-1)));
rep(R1, "R1 = Phi1+Phi2      ");
rep(R2, "R2 = i(Phi1-Phi2)   ");
rep(R3, "R3 = Phi3+Phi4      ");
rep(R4, "R4 = i(Phi3-Phi4)   ");
rep(R1/2+R2, "Phi_D = R1/2+R2     ");
quit;
