default(parisize,"14G");
default(realprecision,220);
M = 900;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
D(s)='q*deriv(s,'q);
ev(s,z)=subst(truncate(s),'q,z);
E6(k)=1-504*sum(m=1,M\k+1,sigma(m,5)*'q^(k*m))+O('q^(M+2));
X = 'q*et(2)*et(16)^2/(et(1)^2*et(8));
T = X/(8*X^2+2*X+1);
PH = -(E6(1)-85*E6(2)+1428*E6(4)-5440*E6(8)+4096*E6(16))/504;
F = PH/D(T);
print("=== level-16 zeta(5) row,  N=16, weight(F)=4, B=2, C=8");
print("x = ", X+O('q^8));
print("t = ", T+O('q^9));
print("Phi_16 = ", PH+O('q^9));
print("F = ", F+O('q^8));
NA=620;
QT = serreverse(T+O('q^(NA+1)));
AT = subst(F+O('q^(NA+1)),'q,QT);
print("A_0..A_10 = ", vector(11,j,polcoef(AT,j-1)));
print("A_n integral for n<=", NA, ": ", denominator(content(truncate(AT)))==1);
S = F/D(X)^2;
print("S = F/(Dx)^2 = ", S+O('q^8));
\\ rational fit of S in X
rf(GG,U,d1,d2,NC)={
  my(cols=List(), Ms);
  for(j=0,d1, listput(cols, Vec(truncate(U^j+O('q^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-GG*U^j+O('q^(NC+1))),-(NC+1))));
  Ms = Mat(vector(#cols,i,cols[i]~));
  matker(Ms);
}
{for(d=0,6, my(k=rf(S,X,d,d,4*d+40)); if(#k>0, print("  S rational in x, degree ",d," kernel: ",k~); break));}
qc = exp(-Pi/2);
sq2 = sqrt(2); l1 = 2+4*sq2; tp = 1/l1;
print("q_c = exp(-2Pi/4) = ", qc);
print("x(q_c) - 1/(2sqrt2) = ", ev(X,qc)-1/(2*sq2));
print("t(q_c) - t_+        = ", ev(T,qc)-tp);
print("Dt(q_c)             = ", ev(D(T),qc));
Fv = ev(F,qc); DFv = ev(D(F),qc); D2T = ev(D(D(T)),qc); Dx = ev(D(X),qc);
print("F(tau_c)            = ", Fv);
print("DF(tau_c)           = ", DFv);
print("DF - (k sqrt N/(4 Pi)) F, k=4 : ", DFv - 4*4/(4*Pi)*Fv);
print("D2t(tau_c)          = ", D2T);
Kf = sqrt(tp*DFv^2/(2*Pi*(-D2T)));
print("K (fold formula)    = ", Kf);
\\ Fricke sign checks
th=1.9; N=16; tau=(cos(th)+I*sin(th))/sqrt(N); taup=-1/(N*tau);
q1=exp(2*Pi*I*tau); q2=exp(2*Pi*I*taup);
print("Phi|_6 W_16 / Phi   = ", ev(PH,q2)/((sqrt(N)*tau)^6*ev(PH,q1)));
print("F|_4 W_16 / F       = ", ev(F,q2)/((sqrt(N)*tau)^4*ev(F,q1)));
print("t(W tau)-t(tau)     = ", ev(T,q2)-ev(T,q1));
print("x(Wtau)*x(tau)      = ", ev(X,q2)*ev(X,q1), "  1/C=",1/8.);
\\ Richardson
MPTS=18; H=25; N0=NA;
Kn(m)=abs(polcoef(AT,m))*tp^m*m^(3/2);
ns=vector(MPTS,j,N0-H*(j-1));
V=matrix(MPTS,MPTS,r,s,1/ns[r]^(s-1)); b=vector(MPTS,r,Kn(ns[r]))~;
sol=matsolve(V,b);
sol2=matsolve(matrix(MPTS-2,MPTS-2,r,s,1/ns[r]^(s-1)),vector(MPTS-2,r,Kn(ns[r]))~);
print("K (Richardson ",MPTS," pts) = ", sol[1]);
print("K (Richardson ",MPTS-2," pts) = ", sol2[1]);
print("K fold / K Richardson - 1  = ", Kf/sol[1]-1);
print("sign of A_n * (-1)^n: ", vector(6,j,sign(polcoef(AT,NA-j+1))));
\\ CM identification, Q(i), tau_c = i/4
W2 = gamma(1/4)^4/Pi^3;
print("Dx(tau_c)           = ", Dx);
print("Dx/W2               = ", Dx/W2, "   algdep 4: ", algdep(Dx/W2,4));
print("(Dx/W2)^2 algdep 4  : ", algdep((Dx/W2)^2,4));
print("K*Pi^(3/2)/W2       = ", Kf*Pi^(3/2)/W2);
print("  algdep deg 2: ", algdep(Kf*Pi^(3/2)/W2,2));
print("  algdep deg 4: ", algdep(Kf*Pi^(3/2)/W2,4));
print("  algdep deg 8: ", algdep(Kf*Pi^(3/2)/W2,8));
print("(K Pi^(3/2)/W2)^2 algdep 4: ", algdep((Kf*Pi^(3/2)/W2)^2,4));
print("K*Pi^(9/2)/Gamma(1/4)^4 = ", Kf*Pi^(9/2)/gamma(1/4)^4);
quit;
