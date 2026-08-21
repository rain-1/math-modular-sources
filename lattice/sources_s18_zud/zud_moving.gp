read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
bpq(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P[n+1],Q[n+1]];
}
N=34; r=zudrow(N); u=r[1]; v=r[2];
print("m : v2(v_m/u_m - v_{m-1}/u_{m-1})   v2(pm(xm)/um - p(m-1)(x(m-1))/u(m-1))   v2((v_m-p_m(x_m))/u_m)");
{ my(pr=0,prev=0,prevz=0);
  for(m=2,N-2,
    my(b=bpq(m,1/2-m), bb=bpq(m-1,1/2-(m-1)));
    my(z=b[1]/u[m+1], zz=bb[1]/u[m]);
    print("  ",m,"   ", valuation(v[m+1]/u[m+1]-v[m]/u[m],2), "        ", valuation(z-zz,2), "        ", valuation((v[m+1]-b[1])/u[m+1],2)));
}
quit;
