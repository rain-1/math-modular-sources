/* 19_abel.gp -- independent check of the cusp periods by Abel summation of the
   Eichler integral D^{-2}Phi along tau = i y (cusp 0, q->1) and tau = 1/2 + i y
   (cusp 1/2, q->-1).  Limits must be the t1- and t2-periods respectively.   */
default(realprecision, 50);
NN = 60000;
re4(n) = my(r=n%5); if(r==1, 1, if(r==4, -1, 0));
im4(n) = my(r=n%5); if(r==2, 1, if(r==3, -1, 0));
cD = vector(NN); cnew3 = vector(NN); cnew4 = vector(NN);
for(n=1, NN, my(X=sumdiv(n,d,re4(d)*d^2)); my(Y=sumdiv(n,d,im4(d)*d^2)); my(U=sumdiv(n,d,re4(n/d)*d^2)); my(V=sumdiv(n,d,im4(n/d)*d^2)); cD[n]=X-2*Y; cnew3[n]=2*U; cnew4[n]=-2*V);
s5=sqrt(5); ph5=(11+5*s5)/2; phm5=(5*s5-11)/2;
cN  = vector(NN, n, cnew3[n] + ph5*cnew4[n]);
cNp = vector(NN, n, cnew3[n] - phm5*cnew4[n]);
S(cv, y, sgn) = my(s=0); my(w=exp(-2*Pi*y)); my(p=1); for(n=1, NN, p=p*w; if(p<1e-45, break); s += cv[n]*sgn^n*p/n^2); s;
z2=zeta(2);
Lch(s) = 5^(-s)*(zetahurwitz(s,1/5)+I*zetahurwitz(s,2/5)-I*zetahurwitz(s,3/5)-zetahurwitz(s,4/5));
LL=Lch(2); LR=real(LL); LI=imag(LL); xi=ph5*LI-LR; xip=-phm5*LI-LR;
ys = [0.02, 0.01, 0.005, 0.0025, 0.00125];
print("cusp 0 (q -> +1), Phi_D  : target zeta(2)/5 = ", z2/5);
for(k=1,5, print("   y=", ys[k], "  S = ", S(cD, ys[k], 1)));
print("cusp 1/2 (q -> -1), Phi_D : target -11*zeta(2)/5 = ", -11*z2/5);
for(k=1,5, print("   y=", ys[k], "  S = ", S(cD, ys[k], -1)));
print("cusp 0 (q -> +1), Phi_new : target xi = ", xi);
for(k=1,5, print("   y=", ys[k], "  S = ", S(cN, ys[k], 1)));
print("cusp 1/2 (q -> -1), Phi'_new : target xi' = ", xip);
for(k=1,5, print("   y=", ys[k], "  S = ", S(cNp, ys[k], -1)));
quit;
