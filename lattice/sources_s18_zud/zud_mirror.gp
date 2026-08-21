read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
NT=44; r=zudrow(NT); u=r[1];
U = vector(NT+1, i, 16^(i-1)*u[i]);   /* integral row */
NN=NT-2;
/* operator L = sum p_i(t) theta^i, from zud_ode2 */
P = [65536*x^4+93184*x^3+384*x^2+28*x, 393216*x^4+491520*x^3+3072*x^2+288*x, 851968*x^4+884736*x^3+23552*x^2+256*x-1, 786432*x^4+655360*x^3+120832*x^2+2688*x+4, 262144*x^4+344064*x^3+112640*x^2+64*x-4];
{ applyL(cv, PP) =
  my(res=vector(NN+1));
  for(i=1,#PP, for(j=0,poldegree(PP[i]),
     my(pc=polcoeff(PP[i],j)); if(pc!=0,
       for(n=0,NN-j, res[n+j+1] += pc*n^(i-1)*cv[n+1]))));
  res;
}
y0 = vector(NN+1,i,U[i]);
print("check L(y0)=0 : ", applyL(y0,P)[1..NN-4]);
/* dL/dtheta */
Pd = vector(4, i, i*P[i+1]);
rhs = -applyL(y0, Pd);   /* need L(h) = rhs */
/* solve for h_n : L(h) coefficient at t^n involves h_n * (indicial at n) ... build triangular solve */
h = vector(NN+1);
{ my(ind = n -> sum(i=1,5, polcoeff(P[i],0)*n^(i-1)));
  for(n=0,NN,
    my(s = 0);
    for(i=1,5, for(j=1,poldegree(P[i]), my(pc=polcoeff(P[i],j));
       if(pc!=0 && n-j>=0, s += pc*(n-j)^(i-1)*h[n-j+1])));
    if(ind(n)==0, h[n+1]=0, h[n+1] = (rhs[n+1]-s)/ind(n)));
}
print("h (first 8): ", vector(8,i,h[i]));
Y0 = sum(n=0,NN,y0[n+1]*x^n)+O(x^(NN+1));
H  = sum(n=0,NN,h[n+1]*x^n)+O(x^(NN+1));
qq = x*exp(H/Y0);
print("q(t) = ", qq+O(x^10));
tq = serreverse(qq);
print("t(q) = ", tq+O(x^16));
print("t(q) coefficients: ", vector(NN-2,i,polcoeff(tq,i)));
print("integral? ", vecmax(vector(NN-2,i,denominator(polcoeff(tq,i)))));
quit;
