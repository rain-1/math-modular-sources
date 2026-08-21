read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
bpq(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P[n+1],Q[n+1]];
}
S(m) = sum(j=1,m,(-1)^j/(2*j-1)^2);
PR=400; z2 = Lp(2,triv,2,PR);
N=34; r=zudrow(N); u=r[1]; v=r[2];
print("m   v2(Psharp/Q - z2)  v2(P/Q - z2)   Psharp-P (exact?)");
{ for(m=1,30,
   my(b=bpq(m,1/2-m));
   my(c = (-1)^m*b[1]/8 - S(m)*u[m+1]);
   print("  m=",m,"  ",valuation(c/u[m+1]-z2,2),"   ",valuation(v[m+1]/u[m+1]-z2,2),"   D=",v[m+1]-c, "  v2(D/Q)=",if(v[m+1]!=c,valuation((v[m+1]-c)/u[m+1],2),"INF")));
}
quit;
