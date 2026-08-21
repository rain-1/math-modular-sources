read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
bpq(n,xx) =
{ my(p0=0,q0=1,p1=1,q1=xx^2-xx+1,p2,q2,a=xx^2-xx);
  if(n==0, return([0,1]));
  for(m=1,n-1,
    q2=((2*m*(m+1)+1+a)*q1-m^2*q0)/(m+1)^2;
    p2=((2*m*(m+1)+1+a)*p1-m^2*p0)/(m+1)^2;
    p0=p1;q0=q1;p1=p2;q1=q2);
  [p1,q1];
}
PR=700; z2=Lp(2,triv,2,PR);
N=70; r=zudrow(N+2); u=r[1]; v=r[2];
G=0;
print("m  v2(Q_m)   -4m+2s2   v2(r_m(x_m))   4m+2-2s2   v2(z2-P/Q)  8m-1-4s2");
{ for(m=1,60,
  G += (-1)^(m-1)/(2*m-1)^2;
  my(b=bpq(m,1/2-m), Tm=(-1)^m*8*(z2-G));
  my(rm = b[1]-Tm*u[m+1], s=hammingweight(m));
  print("  ",m,"  ",valuation(u[m+1],2)," ",-4*m+2*s,"    ",valuation(rm,2)," ",4*m+2-2*s,
        "     ",valuation(v[m+1]/u[m+1]-z2,2)," ",8*m-1-4*s));
}
quit;
