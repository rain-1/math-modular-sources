\\ 05_rows.gp -- (rewritten) the rows on the four-point hosts, one line per pole
\\ placement, including the K-rational placements at a complex special point.
\\ Calibration: Gamma_0(6) pole 1/2 must reproduce Zagier C = CDT's row 1,3,15,93,...
\\ with a THREE-term recurrence of degree 2 and characteristic roots 9, 1.
default(seriesprecision,90);
M = 90;
q = 'q;
i2 = Mod('t,'t^2+1);
w3 = Mod('t,'t^2+3);
E2(n,MM) = 1-24*sum(k=1,MM\n,sigma(k)*q^(n*k))+O(q^MM);
row(F,x) = my(S=subst(F,q,serreverse(x))); vector(M-6,j,polcoeff(S,j-1));
minrec(A) = {my(best=0); for(T=2,6, for(deg=1,6, my(cols=T*(deg+1),rows=40,m,kk); if(rows+T+2>#A,next); m=matrix(rows,cols); for(e=1,rows, my(n=e+1); for(j=0,T-1, for(dd=0,deg, m[e,j*(deg+1)+dd+1]=n^dd*A[n+j+1]))); kk=matker(m); if(#kk>0 && best==0, best=[T,deg,kk[,1]]))); best;}
charroots(b) = {my(T=b[1],deg=b[2],v=b[3],p=0); p=sum(j=0,T-1, lift(v[(j+1)*(deg+1)])*'y^j); p;}
show(name,A) = {my(b=minrec(A)); print(name); print("   A_0..A_5 = ",vector(6,i,A[i])); if(b==0, print("   NO recurrence with <=6 terms, deg<=6"), print("   minimal recurrence: ",b[1]," terms, degree ",b[2],"   CDT-shape(3 terms)? ",if(b[1]==3,"YES","NO")); print("   characteristic polynomial: ",charroots(b)));}

print("################ CALIBRATION: Gamma_0(6), weight 1, chi_-3 ################");
G6 = znstar(6,1);
mf6 = mfinit([6,1,[G6,[1]]],3);
Bs = mfbasis(mf6);
F6 = 3*(Ser(mfcoefs(Bs[1],M),q)+Ser(mfcoefs(Bs[2],M),q));
h6 = q*eta(q^2+O(q^M))*eta(q^6+O(q^M))^5/(eta(q+O(q^M))^5*eta(q^3+O(q^M)));
show("Gamma_0(6) pole cusp 1/2 = Zagier C = CDT  (expect 1,3,15,93; 3 terms, deg 2; roots 9,1)", row(F6,h6/(1+9*h6)));
show("Gamma_0(6) pole cusp 1/3 = Zagier A", row(F6,h6/(1+8*h6)));
show("Gamma_0(6) pole cusp 0   = Zagier F", row(F6,h6));

print();
print("################ X_0(5): the Q(i)-rational placement ################");
h5 = q*eta(q^5+O(q^M))^6/eta(q+O(q^M))^6;
F5 = (5*E2(5,M)-E2(1,M))/4;
print("special points: mu(cusp 0)=0, mu(ell_2 pair) = -11 +- 2i  (02_geom.out)");
print("pole at cusp 0     : lambda = -11 +- 2i, |lambda_1|=|lambda_2|=sqrt(125): DEAD (no dominant root)");
print("pole at ell_2      : lambda = 11-2i (|.|=11.1803) and -4i (|.|=4), K=Q(i), N(lambda_2)=16 (not a unit)");
show("  X_0(5) weight-2 row, pole at ell_2, K=Q(i)", row(F5,h5/(1-(-11+2*i2)*h5)));
show("  X_0(5) weight-2 row, pole at cusp 0", row(F5,h5));

print();
print("################ X_0(7): the Q(sqrt-3)-rational placement ################");
h7 = q*eta(q^7+O(q^M))^4/eta(q+O(q^M))^4;
F7e = (7*E2(7,M)-E2(1,M))/6;
G7 = znstar(7,1);
mf7 = mfinit([7,1,[G7,[3]]],3);
F7 = Ser(mfcoefs(mfbasis(mf7)[1],M),q);
print("special points: mu(cusp 0)=0, mu(ell_3 pair) = -13/2 +- (3 sqrt3/2) i");
print("pole at cusp 0 : lambda = -13/2 +- (3sqrt3/2)i, both |.|=7: DEAD");
print("pole at ell_3  : lambda = (13-3sqrt-3)/2 (|.|=7) and -3sqrt-3 (|.|=5.196), K=Q(sqrt-3), N(lambda_2)=27");
show("  X_0(7) weight-1 row (nebentypus chi_-7), pole at ell_3, K=Q(sqrt-3)", row(F7,h7/(1-((-13+3*w3)/2)*h7)));
show("  X_0(7) weight-1 row (nebentypus chi_-7), pole at cusp 0", row(F7,h7));
show("  X_0(7) weight-2 row (trivial nebentypus), pole at ell_3, K=Q(sqrt-3)", row(F7e,h7/(1-((-13+3*w3)/2)*h7)));

print();
print("################ X_0(9): both conjugate placements ################");
h9 = q*eta(q^9+O(q^M))^3/eta(q+O(q^M))^3;
G9 = znstar(9,1);
mf9 = mfinit([9,1,[G9,[3]]],3);
Bs9 = mfbasis(mf9);
F9 = Ser(mfcoefs(Bs9[1],M),q);
print("special points: mu(cusp 0)=0, mu(cusp 1/3), mu(cusp 2/3) = -9/2 -+ (3sqrt3/2) i");
print("pole at cusp 0   : lambda = -9/2 +- (3sqrt3/2)i, |.|=sqrt27 both: DEAD");
print("pole at cusp 1/3 : lambda = (9+3sqrt-3)/2 (|.|=sqrt27) and 3sqrt-3 (|.|=sqrt27): DEAD, equal moduli");
show("  X_0(9) weight-1 row, pole at cusp 0 = Zagier B", row(F9,h9));
show("  X_0(9) weight-1 row, pole at cusp 1/3, K=Q(sqrt-3)", row(F9,h9/(1-((-9-3*w3)/2)*h9)));
quit;
