/* 03_ident.gp -- TASK 3.  Identification: F and F^2 in M_*(Gamma_0(N)); Phi is NOT
   holomorphic; Atkin-Lehner signs; the critical L-value L(Xi,2).                */
default(parisize, 1000000000);
read("lib.gp");
P = 200;
{ e4(d,P) = my(s=1+O('q^P)); for(n=1,(P-1)\d, s += 240*sigma(n,3)*'q^(d*n)); s; }
{ etap(d,e,P) = my(t=1+O('q^P)); for(n=1,P\d, t *= (1-'q^(d*n))); 'q^(d*e/24)*t^e; }
{ nf(Nl,P) = my(mf=mfinit([Nl,4],0), fs=mfeigenbasis(mf)); Ser(mfcoefs(fs[1],P-1),'q); }
{ ident(f,Nl,wt,nfit,ncheck,lab) = my(mf,vv,co,mfb,chk);
  mf = mfinit([Nl,wt],4); vv = vector(nfit,j,polcoeff(f,j-1));
  co = mftobasis(mf,vv,1); mfb = mfcoefs(mflinear(mf,co),ncheck);
  chk = vector(ncheck+1,j, mfb[j]-polcoeff(f,j-1));
  print(lab,"  dim=",mfdim(mf),"  coords=",co,"   max defect over ",ncheck," coeffs: ",vecmax(apply(abs,chk))); }
{ for(k=1,3, my(R=ROWS[k], S=Setup(k,P));
  ident(S[2],  R[2],2,40,120,concat(R[1]," F   in M_2(Gamma_0(N)) "));
  ident(S[2]^2,R[2],4,60,150,concat(R[1]," F^2 in M_4(Gamma_0(N)) "));
  ident(S[4],  R[2],4,60,150,concat(R[1]," Phi in M_4(Gamma_0(N)) "))); }
print();
print("=== explicit decompositions of F^2 (exact to O(q^199)) ===");
S1=Setup(1,P); S2=Setup(2,P); S3=Setup(3,P);
f7=nf(7,P); f5=nf(5,P); f5b=subst(f5,'q,'q^2)+O('q^P);
g9=etap(3,8,P); g9b=etap(6,8,P);
print(" s7 : F^2 - [(E4+49E4(7t))/50 + (16/5) f_{7.4.a.a}] = ", S1[2]^2-((e4(1,P)+49*e4(7,P))/50+16/5*f7));
print(" s10: F^2 - [(E4+4E4(2t)+25E4(5t)+100E4(10t))/130 + (28/13)(f_{5.4.a.a}+4 f_{5.4.a.a}(2t))] = ", S2[2]^2-((e4(1,P)+4*e4(2,P)+25*e4(5,P)+100*e4(10,P))/130+28/13*(f5+4*f5b)));
print(" s18: F^2 - [(E4(3t)+4E4(6t))/5 + 12(eta(3t)^8+4 eta(6t)^8)] = ", S3[2]^2-((e4(3,P)+4*e4(6,P))/5+12*(g9+4*g9b)));
print();
\p 50
NT = 400;
{ E2n(t) = my(q=exp(2*Pi*I*t),s=1,p=1); for(n=1,NT, p*=q; s -= 24*sigma(n)*p); s; }
u18(t) = (eta(2*t,1)*eta(3*t,1)^2*eta(18*t,1)/(eta(t,1)*eta(6*t,1)^2*eta(9*t,1)))^6;
F18(t) = (-6*E2n(t)+12*E2n(2*t)+36*E2n(3*t)-72*E2n(6*t)-54*E2n(9*t)+108*E2n(18*t))/24;
x18(t) = my(u=u18(t)); u/(1+14*u+u^2);
Ph18(t) = my(u=u18(t),F=F18(t)); F^2*u*(1-u^2)/(1+14*u+u^2)^2;
W9(t) = (9*t-5)/(18*t-9);
tt = 0.211+0.4013*I;
print("=== s18 Atkin-Lehner W_9 (fixed point tau_0=(3+i)/6) ===");
print("  u|W_9 / u        = ", u18(W9(tt))/u18(tt));
print("  x|W_9 / x        = ", x18(W9(tt))/x18(tt));
print("  F|_2 W_9 / F     = ", (F18(W9(tt))/(6*tt-3)^2)/F18(tt), "   [WEIGHT_DROP 4.3 asserts -1]");
print("  Phi|_4 W_9 / Phi = ", (Ph18(W9(tt))/(6*tt-3)^4)/Ph18(tt));
print("  u|W_18 * u       = ", u18(-1/(18*tt))*u18(tt));
print("  F|_2 W_18 / F    = ", (F18(-1/(18*tt))/(sqrt(18)*tt)^2)/F18(tt));
print();
\p 60
M=400; NN=M+6;
{ CV=vector(3); for(k=1,3, my(S=Setup(k,NN)); CV[k]=vector(M,m,polcoeff(S[4],m))); }
{ Ifun(k,s,Nl)=my(t=0); for(m=1,M, t += CV[k][m]*(2*Pi*m)^(-s)*incgam(s,2*Pi*m/sqrt(Nl))); t; }
{ Lam(k,s,Nl)=Ifun(k,s,Nl)-Nl^(2-s)*Ifun(k,4-s,Nl); }
print("=== critical values (Xi = D^-1 Phi = sum c'(m) q^m) ===");
{ for(k=1,3, my(Nl=ROWS[k][2], xi);
  xi = 4*Pi^3*Lam(k,3,Nl);
  print("  ",ROWS[k][1],"  Lambda(Phi,2) = ",Lam(k,2,Nl));
  print("        L(Xi,2) = ",xi,"   / zeta(2) = ",xi/zeta(2),"   / L(chi_-3,2) = ",xi/lfun(-3,2))); }
quit;
