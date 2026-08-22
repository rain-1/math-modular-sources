/* J_m(s) = int_0^1 int_0^1 [x(1-x)y(1-y)]^m (xy)^{s-1} /(1-xy)^{m+1} dx dy
   = C_m(s)*zeta(2,s) - Dm(s), C_m,Dm in Q.
   R(k) = m! prod_{j=1}^m (k+j) / prod_{i=m}^{2m}(k+s+i)^2 . */
Jform(m,s)=
{
  my(R, A=vector(m+1), B=vector(m+1), num, CA, CD, poles);
  poles = vector(m+1, t, m+t-1);       /* i = m..2m */
  num(k) = m! * prod(j=1,m, k+j);
  /* A_i = num(-(s+i)) / prod_{i'!=i} (i'-i)^2  ; B_i = d/dk [ num(k)/prod_{i'!=i}(k+s+i')^2 ] at k=-(s+i) */
  for(t=1,m+1,
    my(i=poles[t], k0=-(s+i), pr, dpr, f, df);
    pr = prod(t2=1,m+1, if(t2==t,1,(poles[t2]-i)^2));
    f = num(k0)/pr;
    A[t] = f;
    /* logarithmic derivative: d/dk log(num) - 2*sum_{i'!=i} 1/(k+s+i') at k0 */
    df = sum(j=1,m, 1/(k0+j)) - 2*sum(t2=1,m+1, if(t2==t,0,1/(poles[t2]-i)));
    B[t] = f*df;
  );
  CA = sum(t=1,m+1, A[t]);
  if(abs(sum(t=1,m+1,B[t]))>0, error("sum B nonzero"));
  CD = sum(t=1,m+1, A[t]*sum(l=0,poles[t]-1, 1/(l+s)^2)) + sum(t=1,m+1, B[t]*sum(l=0,poles[t]-1, 1/(l+s)));
  [CA, CD];
}
{
for(m=0,6,
  my(u=Jform(m,1/3), v=Jform(m,2/3));
  print("m=",m,"  C(1/3)=",u[1],"  C(2/3)=",v[1],"  ratio=",u[1]/v[1]);
);
}
