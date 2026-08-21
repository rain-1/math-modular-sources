read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
/* Beukers Theta-Pade: (n+1)^2 w_{n+1} = (2n(n+1)+1-x+x^2) w_n - n^2 w_{n-1} */
bpq(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P[n+1],Q[n+1]];
}
N=14; r=zudrow(N); u=r[1]; v=r[2];
print("check q_n(x) at fixed x=-n+1/2 (moving): ");
{ for(n=0,10, my(b=bpq(n,-n+1/2)); print("  n=",n,"  q_n=",b[2],"   u_n=",u[n+1],"   equal? ",b[2]==u[n+1])); }
print("p_n(-n+1/2) vs v_n:");
{ for(n=0,8, my(b=bpq(n,-n+1/2)); print("  n=",n,"  p_n=",b[1],"  v_n=",v[n+1],"  v_n-p_n=",v[n+1]-b[1],"  (v_n-p_n)/u_n=",if(u[n+1]!=0,(v[n+1]-b[1])/u[n+1]))); }
print("Beukers at fixed x=1/2 (Calegari's approximations): p_n/q_n, 2-adic");
{ my(z2); for(n=2,12,my(b=bpq(n,1/2)); print("  n=",n," p_n/q_n = ",b[1]/b[2]," v_2 num/den ",valuation(b[1],2),",",valuation(b[2],2))); }
quit;
