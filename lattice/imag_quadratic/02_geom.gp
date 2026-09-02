\\ 02_geom.gp -- exact special-point configuration of every four-point genus-zero
\\ Gamma_H host (N<=60 list from 01_groups.py: X_0(5),X_0(6),X_0(7),X_0(8),X_0(9),X_1(5)).
\\ Reference Hauptmodul h with h = q + O(q^2) at the cusp infinity; mu := 1/h.
\\ Pole placement at special point P_j gives lambda_i = mu_i - mu_j  (hostscan REPORT sec.3).
default(realprecision,80);
Y = 60;
cuspt(a,b,c,d,Y) = (a*I*Y+b)/(c*I*Y+d);
h5(t) = (eta(5*t,1)/eta(t,1))^6;
h6(t) = eta(2*t,1)*eta(6*t,1)^5/(eta(t,1)^5*eta(3*t,1));
h7(t) = (eta(7*t,1)/eta(t,1))^4;
h8(t) = eta(4*t,1)^5*eta(8*t,1)^2/(eta(t,1)^2*eta(2*t,1)^5);
h9(t) = (eta(9*t,1)/eta(t,1))^3;
rec(z) = [bestappr(real(z),10^20), bestappr(imag(z)^2,10^20)];
show(nm, z) = print(nm, " = ", z, "   [Re ~ ", bestappr(real(z),10^15), " , (Im)^2 ~ ", bestappr(imag(z)^2,10^15), "]");

print("=== X_0(5): cusps oo,0 ; two order-2 elliptic points ===");
print("h5 at cusp 0 -> infinity, so mu(cusp 0)=0 (checked: |h5(i/60)| large)");
te1 = (-2+I)/5;
te2 = (2+I)/5;
show("mu(ell_1)", 1/h5(te1));
show("mu(ell_2)", 1/h5(te2));

print();
print("=== X_0(6): cusps oo,0,1/2,1/3 ===");
show("mu(cusp 0)  ", 1/h6(cuspt(0,-1,1,0,Y)));
show("mu(cusp 1/2)", 1/h6(cuspt(1,0,2,1,Y)));
show("mu(cusp 1/3)", 1/h6(cuspt(1,0,3,1,Y)));

print();
print("=== X_0(7): cusps oo,0 ; two order-3 elliptic points ===");
show("mu(cusp 0)  ", 1/h7(cuspt(0,-1,1,0,Y)));
show("mu(ell_1)", 1/h7((-5+I*sqrt(3))/14));
show("mu(ell_2)", 1/h7((-9+I*sqrt(3))/14));

print();
print("=== X_0(8): cusps oo,0,1/2,1/4 ===");
show("mu(cusp 0)  ", 1/h8(cuspt(0,-1,1,0,Y)));
show("mu(cusp 1/2)", 1/h8(cuspt(1,0,2,1,Y)));
show("mu(cusp 1/4)", 1/h8(cuspt(1,0,4,1,Y)));

print();
print("=== X_0(9): cusps oo,0,1/3,2/3 ===");
show("mu(cusp 0)  ", 1/h9(cuspt(0,-1,1,0,Y)));
show("mu(cusp 1/3)", 1/h9(cuspt(1,0,3,1,Y)));
show("mu(cusp 2/3)", 1/h9(cuspt(2,1,3,2,Y)));
quit;
