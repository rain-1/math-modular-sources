/* 03_model.gp -- Weierstrass model of the K3 row from the J-map,
 * and the check that the Apery holomorphic solution is the period.
 */
default(parisizemax, 4000000000);
NT = 40;

uu = -644972544*x^15 - 89579520*x^14 + 569161728*x^13 - 161970688*x^12 - 152715264*x^11 + 120155136*x^10 - 22576640*x^9 - 11678400*x^8 + 9492480*x^7 - 3411456*x^6 + 787680*x^5 - 126072*x^4 + 14080*x^3 - 1056*x^2 + 48*x - 1;
vv = x^4*(x-1)^6*(27*x^2-10*x+1)^3;
p4 = 216*x^4 + 64*x^3 - 48*x^2 + 12*x - 1;
q8 = 5832*x^8 + 34560*x^7 - 30016*x^6 + 12624*x^5 - 4380*x^4 + 1280*x^3 - 240*x^2 + 24*x - 1;
print("U == -(4t-1)^3 P4^3 : ", uu == -(4*x-1)^3*p4^3);
print("U - 1728 V == -Q8^2 : ", uu - 1728*vv == -q8^2);

c4p = (4*x-1)*p4;
c6p = -q8;
print("c4^3 - c6^2 == -1728 V : ", c4p^3 - c6p^2 == -1728*vv);
print("c4 = ", c4p);
print("c6 = ", c6p);
print("Delta = ", -vv);

/* --- the Apery series and the canonical nome --------------------------- */
uvec = vector(NT+2); gvec = vector(NT+2);
uvec[1] = 1; gvec[1] = 0;
aa=11; bb=11; cc=4; dd=37; ee=0; ff=3; gg=27; hh=-27; jj=6;
/* P(n)=11n^2+11n+4, Q(n)=37n^2+3, R(n)=3(3n-1)(3n-2)=27n^2-27n+6 */
{
for(n=1,NT+1,
  my(n1=n-1, P0, Q0, R0, dP, dQ, dR);
  P0 = aa*n1^2+bb*n1+cc; Q0 = dd*n1^2+ee*n1+ff; R0 = gg*n1^2+hh*n1+jj;
  dP = 2*aa*n1+bb; dQ = 2*dd*n1+ee; dR = 2*gg*n1+hh;
  uvec[n+1] = (P0*uvec[n] - if(n>=2,Q0*uvec[n-1],0) + if(n>=3,R0*uvec[n-2],0))/n^2;
  gvec[n+1] = (dP*uvec[n]+P0*gvec[n] - if(n>=2, dQ*uvec[n-1]+Q0*gvec[n-1],0) + if(n>=3, dR*uvec[n-2]+R0*gvec[n-2],0) - 2*n*uvec[n+1])/n^2;
);
}
print("u_n = ", vector(12,i,uvec[i]));
y0 = sum(n=0,NT, uvec[n+1]*x^n) + O(x^(NT+1));
gs = sum(n=0,NT, gvec[n+1]*x^n) + O(x^(NT+1));
qq = x*exp(gs/y0);
print("nome q = ", qq + O(x^8));
QQ = -qq^4;

e4 = 1 + 240*sum(k=1,NT, sigma(k,3)*z^k) + O(z^(NT+1));
e6 = 1 - 504*sum(k=1,NT, sigma(k,5)*z^k) + O(z^(NT+1));
E4 = subst(e4, z, QQ);
E6 = subst(e6, z, QQ);

print();
print("E4(Q)/y0^4 (first terms) = ", E4/y0^4 + O(x^8));
print("E4(Q)*y0^4 (first terms) = ", E4*y0^4 + O(x^8));
print("target c4 = ", c4p);
print();
r1 = (E4/y0^4)/c4p; print("ratio (E4/y0^4)/c4 = ", r1 + O(x^8));
r2 = (E4*y0^4)/c4p; print("ratio (E4*y0^4)/c4 = ", r2 + O(x^8));
print();
s1 = (E6/y0^6)/c6p; print("ratio (E6/y0^6)/c6 = ", s1 + O(x^8));
s2 = (E6*y0^6)/c6p; print("ratio (E6*y0^6)/c6 = ", s2 + O(x^8));
quit;
