/* Exact evaluation of the cyclotomic cellular integrals of SOL note 2, section 5:
   C_n = int_0^1 int_0^1 x^{2n}(1-x^4)^n y^{2n}(1-y^4)^n /(1+x^2y^2)^{2n+1} dx dy.
   Period basis: [q0,q1,q2]  <->  q0 + q1*Pi + q2*G.
   Reduction: u=xy gives
     int int x^A y^B/(1+x^2y^2)^r = (Jv(B,r)-Jv(A,r))/(A-B)   (A != B)
                                  = Kv(A,r)                    (A == B)
   with Jv(m,r)=int_0^1 u^m/(1+u^2)^r du, Kv(m,r)=-int_0^1 u^m log u/(1+u^2)^r du.
   Avoids PARI builtin names psi,M,Phi,S,cmp. */

Jz = List(); Kz = List();   /* Jz[r+1] = Jv(0,r), Kz[r+1] = Kv(0,r) */
initzero(rmax) = {
  Jz = vector(rmax+1); Kz = vector(rmax+1);
  Jz[1] = [1,0,0];            /* Jv(0,0) = 1 */
  Kz[1] = [1,0,0];            /* Kv(0,0) = 1 */
  Jz[2] = [0,1/4,0];          /* pi/4 */
  Kz[2] = [0,0,1];            /* G */
  for(r=2,rmax,
    Jz[r+1] = [1/(2^r*(r-1)),0,0] + (2*r-3)/(2*r-2)*Jz[r];
    Kz[r+1] = (Jz[r] + (2*r-3)*Kz[r])/(2*r-2);
  );
}

memJ = 0; memK = 0;
initmem(mmax,rmax) = {
  memJ = matrix(mmax\2+1, rmax+1); memK = matrix(mmax\2+1, rmax+1);
  for(r=0,rmax, memJ[1,r+1] = Jz[r+1]; memK[1,r+1] = Kz[r+1]);
  for(s=1,mmax\2, for(r=0,rmax,
    m = 2*s;
    if(r==0,
      memJ[s+1,1] = [1/(m+1),0,0]; memK[s+1,1] = [1/(m+1)^2,0,0],
      memJ[s+1,r+1] = memJ[s,r] - memJ[s,r+1];
      memK[s+1,r+1] = memK[s,r] - memK[s,r+1];
    );
  ));
}
Jv(m,r) = memJ[m\2+1, r+1];
Kv(m,r) = memK[m\2+1, r+1];

cellular(n) = {
  my(r=2*n+1, tot=[0,0,0], c, A, B);
  for(i=0,n, for(j=0,n,
    c = (-1)^(i+j)*binomial(n,i)*binomial(n,j);
    A = 2*n+4*i; B = 2*n+4*j;
    if(A==B, tot += c*Kv(A,r), tot += c*(Jv(B,r)-Jv(A,r))/(A-B));
  ));
  tot;
}

NMAX = eval(getenv("NMAX")); if(NMAX==0, NMAX=14);
initzero(2*NMAX+2); initmem(6*NMAX+4, 2*NMAX+2);
Uv = vector(NMAX+1); Vv = vector(NMAX+1);
{
print("n | pi-coeff (must be 0) | A_n | B_n | U_n | V_n");
for(nn=0,NMAX,
  my(v=cellular(nn), An=v[3], Bn=v[1], pic=v[2], u, w);
  u = (-1)^nn*16^nn*An; w = (-1)^nn*16^nn*Bn;
  Uv[nn+1]=u; Vv[nn+1]=w;
  print(nn," pi=",pic,"  U=",u,"  V=",w,
        "  denU=",denominator(u),"  denV=",denominator(w),
        "  -V/U-G=",if(u!=0, -w/u-Catalan, 0));
);
}
write("out_cellular.txt", "Uv=",Uv);
write("out_cellular.txt", "Vv=",Vv);
quit();
