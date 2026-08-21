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
N=20; r=zudrow(N); u=r[1]; v=r[2];
print("m  P_m   conj = (-1)^(m+1)*p_m(x_m)/8 - S_m*Q_m   ratio/diff");
{ for(m=0,16,
   my(b=bpq(m,1/2-m));
   my(c = (-1)^(m+1)*b[1]/8 - S(m)*u[m+1]);
   print("  m=",m,"  P=",v[m+1],"  conj=",c,"  diff=",v[m+1]-c, "  P/conj=",if(c!=0,v[m+1]/c)));
}
quit;
