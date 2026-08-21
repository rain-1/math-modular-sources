\\ identify eta-quotient parametrisation of T=(12,4,16):  t Hauptmodul, F=sum a_n t^n
read("/home/ubuntu/code/math-modular-sources/lattice/zeta3_lattice/rows.gp");
T=rowsI(12,4,16); A=T[1];
PR=18;
E(d)=my(s=1+O(q^PR));for(n=1,(PR-1)\d,s*=(1-q^(d*n)));s;   \\ eta(d tau) q^{-d/24}
LE=[E(1),E(2),E(4),E(8)]; LL=vector(4,i,log(LE[i]));
{for(r4=-8,8,for(r8=-8,8,
  my(r2=24-3*r4-7*r8, r1=-(r2+r4+r8), r=[r1,r2,r4,r8]);
  if(vecmax(abs(r))>16,next);
  my(t=q*prod(i=1,4,LE[i]^r[i])+O(q^PR));
  my(F=0*q+O(q^PR),tp=1+O(q^PR));
  for(n=0,PR-1, F+=A[n+1]*tp; tp*=t);
  my(lf=log(F+O(q^PR)));
  my(M=matrix(4,4,i,j,polcoeff(LL[j],i)), v=vector(4,i,polcoeff(lf,i))~);
  if(matdet(M)==0,next);
  my(s=matsolve(M,v));
  if(denominator(s)!=1,next);
  my(chk=lf-sum(i=1,4,s[i]*LL[i]));
  if(truncate(chk)==0,
    print("t = q * prod eta(d tau)^r_d /q^{...}, r(1,2,4,8)=",r,
          "   F = prod eta(d tau)^s_d, s=",Vec(s~)," weight=",vecsum(Vec(s~))/2))));}
\q
