read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
bpq(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P[n+1],Q[n+1]];
}
S(m) = sum(j=1,m,(-1)^j/(2*j-1)^2);
N=34; r=zudrow(N); u=r[1]; v=r[2];
print("m  v2(Q_m) v2(P_m)  v2(pm) D_m=P_m-conj  v2(D_m)  v2(D_m/Q_m)  |D_m/Q_m|");
{ for(m=1,30,
   my(b=bpq(m,1/2-m));
   my(c = (-1)^(m+1)*b[1]/8 - S(m)*u[m+1]);
   my(D = v[m+1]-c);
   print("  m=",m,"  vQ=",valuation(u[m+1],2)," vP=",valuation(v[m+1],2)," vp=",valuation(b[1],2)," vD=",valuation(D,2)," vD/Q=",valuation(D/u[m+1],2),"  |D/Q|=",abs(D/u[m+1])*1.0));
}
quit;
