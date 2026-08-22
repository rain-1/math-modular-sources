/* mirror map / nome for a MUM row  L = P th^2 + Q th + R,
   P = 1 - al t + de t^2,  Q = -be t + q2 t^2,  R = -ga t + r2 t^2.
   y1 = sum a_n t^n ;  y2 = y1 log t + h ;  L h = -(2 P th(y1) + Q y1) ;
   q = t exp(h/y1) ; t(q) = serreverse(q).                                   */
mirror(al,be,ga,de,q2,r2,M) =
{
  my(a=vector(M+3)); a[1]=1; a[2]=ga;
  for(n=1,M+1, my(P=al*n^2+be*n+ga, Q=de*n^2+(-2*de+q2)*n+(de-q2+r2));
    a[n+2]=(P*a[n+1]-Q*a[n])/(n+1)^2);
  my(y1=sum(n=0,M,a[n+1]*'t^n)+O('t^(M+1)));
  my(Pp=1-al*'t+de*'t^2, Qp=-be*'t+q2*'t^2, Rp=-ga*'t+r2*'t^2);
  my(src = -(2*Pp*('t*deriv(y1,'t)) + Qp*y1));
  my(h=vector(M+2));   \\ h[i] = h_{i-1}, h_0 = 0
  for(n=1,M,
    my(s=polcoeff(src,n));
    my(c1 = if(n>=1, (al*(n-1)^2 + be*(n-1) + ga)*h[n], 0));
    my(c2 = if(n>=2, (de*(n-2)^2 + q2*(n-2) + r2)*h[n-1], 0));
    h[n+1] = (s + c1 - c2)/n^2);
  my(H=sum(n=1,M,h[n+1]*'t^n)+O('t^(M+1)));
  my(res = Pp*('t*deriv('t*deriv(H,'t),'t)) + Qp*('t*deriv(H,'t)) + Rp*H - src);
  print("  L h - src residual (should be O(t^M)): ", if(res==0, "0", Ser(vector(4,i,polcoeff(res,i-1)),'t)));
  my(qq='t*exp(H/y1));
  [qq, serreverse(qq), y1];
};
{
default(realprecision,60);
print("=== sanity: Apery sqrt row (136,68,10,16,-16,4) -> q2=de=16, r2=ze=4 ===");
r=mirror(136,68,10,16,16,4,20);
print("  t(q) = ", Ser(vector(9,i,polcoeff(r[2],i-1)),'q));
print("=== new row (234,-39,-78,-27,63,-30):  q2 = ep+2de = 9, r2 = ze-de+q2 = 6 ===");
r=mirror(234,-39,-78,-27,9,6,32);
print("  q(t) = ", Ser(vector(8,i,polcoeff(r[1],i-1)),'t));
print("  t(q) = ", Ser(vector(12,i,polcoeff(r[2],i-1)),'q));
print("  t(q) integral? ", vecmax(vector(30,i,denominator(polcoeff(r[2],i))))==1);
print("  denominators: ", vector(14,i,denominator(polcoeff(r[2],i))));
}
