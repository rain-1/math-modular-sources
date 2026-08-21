read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
bpqall(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
N=20; r=zudrow(N); u=r[1]; v=r[2];
b=bpqall(N,1/2);
{for(n=0,14, print("n=",n,"  q_n(1/2)=",b[2][n+1],"  Q_n=",u[n+1],"  q/Q=",b[2][n+1]/u[n+1],
   "   p_n(1/2)/8 =",b[1][n+1]/8, "  P_n=",v[n+1], "  ratio ",if(v[n+1]!=0,b[1][n+1]/8/v[n+1]),
   "   diff of ratios p/(8q) - P/Q = ", b[1][n+1]/(8*b[2][n+1]) - v[n+1]/u[n+1]));}
quit;
