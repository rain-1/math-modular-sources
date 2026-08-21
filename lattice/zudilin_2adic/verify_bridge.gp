/* Exact verification of the bridge identity to large m.
   Q_m = q_m(1/2-m),  P_m = (-1)^m p_m(1/2-m)/8 + G_m Q_m,  G_m = sum_{j<=m} (-1)^(j-1)/(2j-1)^2 */
read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
default(parisizemax, 6000000000);
bpq(n,xx) =
{ my(p0=0,q0=1,p1=1,q1=xx^2-xx+1,p2,q2,a=xx^2-xx);
  if(n==0, return([0,1])); 
  for(m=1,n-1,
    q2=((2*m*(m+1)+1+a)*q1-m^2*q0)/(m+1)^2;
    p2=((2*m*(m+1)+1+a)*p1-m^2*p0)/(m+1)^2;
    p0=p1;q0=q1;p1=p2;q1=q2);
  [p1,q1];
}
M = 160;
r=zudrow(M+2); u=r[1]; v=r[2];
G=0; okq=1; okp=1; lastq=0; lastp=0;
{ for(m=1,M,
   G += (-1)^(m-1)/(2*m-1)^2;
   my(b=bpq(m,1/2-m));
   if(b[2]!=u[m+1], okq=0; if(lastq==0,lastq=m));
   if((-1)^m*b[1]/8 + G*b[2] != v[m+1], okp=0; if(lastp==0,lastp=m));
  );
}
print("Q_m = q_m(1/2-m) for 1<=m<=",M,": ", okq, "  first failure ", lastq);
print("P_m = (-1)^m p_m(1/2-m)/8 + G_m Q_m for 1<=m<=",M,": ", okp, "  first failure ", lastp);
quit;
