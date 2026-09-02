default(realprecision,70);
K=0.68537184849053440595143004488054241747925748918627;
M=260;
et(k)=prod(n=1,M+1,1-q^(k*n)+O(q^(M+2)));
h=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/q;
X=h/((h+3)*(h+4));
E4(k)=1+240*sum(n=1,M+1,sigma(n,3)*q^(k*n))+O(q^(M+2));
E=(-9*E4(3)+16*E4(4))/7;
D(s)=q*deriv(s,q);
qc=exp(-Pi/sqrt(3));   \\ q at tau_* = i/sqrt(12)
ev(s)=subst(truncate(s),q,qc);
print("x(qc) = ",ev(X),"   x_+ = ",7-4*sqrt(3));
print("Dx(qc) = ",ev(D(X)),"  (should be ~0)");
D2x=ev(D(D(X)));  Ev=ev(E); DEv=ev(D(E)); hc=ev(h);
print("h_12(qc) = ",hc,"  (2sqrt3=",2*sqrt(3),")");
Qc=143*(7-4*sqrt(3))^3+189*(7-4*sqrt(3))^2+21*(7-4*sqrt(3))+7;
\\ G = 7E/Q(x); G'(q) = 7 DE/(q Q) at the fold (x'=0). x''(qc) = D2x/qc^2.
Gp=7*DEv/(qc*Qc);
xpp=D2x/qc^2;
hx=Gp*sqrt(2*(7-4*sqrt(3))/(-xpp));
Kf=hx/(2*sqrt(Pi));
print("K (Richardson) = ",K);
print("K (fold formula)= ",Kf, "   ratio ",Kf/K);
print("K^2 = ",K^2, "   49 x+ (DE)^2/(2 pi Q^2 (-D2x)) = ", 49*(7-4*sqrt(3))*DEv^2/(2*Pi*Qc^2*(-D2x)));
\\ CM values for later identification
print("E(tau*) = ",Ev,"  DE(tau*) = ",DEv,"  D2x(tau*) = ",D2x);
E2(k)=1-24*sum(n=1,M+1,sigma(n,1)*q^(k*n))+O(q^(M+2)); E6(k)=1-504*sum(n=1,M+1,sigma(n,5)*q^(k*n))+O(q^(M+2));
print("E4(3t*),E4(4t*),E6(3t*),E6(4t*),E2(3t*),E2(4t*),E2(t*): ",[ev(E4(3)),ev(E4(4)),ev(E6(3)),ev(E6(4)),ev(E2(3)),ev(E2(4)),ev(E2(1))]);
write("cmvals.txt",[K,Ev,DEv,D2x,hc,ev(E4(3)),ev(E4(4)),ev(E6(3)),ev(E6(4)),ev(E2(3)),ev(E2(4)),ev(E2(1))]);
quit;
