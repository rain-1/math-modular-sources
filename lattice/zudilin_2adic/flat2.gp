read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
bpqall(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
bpq(n,xx)=my(b=bpqall(n,xx));[b[1][n+1],b[2][n+1]];
S(m) = sum(j=1,m,(-1)^j/(2*j-1)^2);
PR=400; z2 = Lp(2,triv,2,PR);
/* verify shift equation numerically with high-n approximants */
NN=45;
print("shift check: Theta(x)+Theta(x+1)+2/x^2 valuations");
{ for(m=1,5, my(a=bpq(NN,1/2-m), c=bpq(NN,1/2-m+1), xx=1/2-m);
    my(th1=a[1]/a[2], th0=c[1]/c[2]);
    print("  m=",m,"  v2(Th(x_m)+Th(x_{m-1})+2/x_m^2) = ", valuation(th1+th0+2/xx^2,2)));}
print("");
print("Theta_m vs (-1)^m*8*(z2 - S_m):");
{ for(m=0,6, my(a=bpq(NN,1/2-m)); print("  m=",m," v2(diff)=",valuation(a[1]/a[2] - (-1)^m*8*(z2-S(m)),2)));}
print("");
N=34; r=zudrow(N); u=r[1]; v=r[2];
print("m   v2(Pflat/Q - z2)   v2(P/Q - z2)   v2(P/Q - Pflat/Q)");
{ for(m=1,30,
   my(b=bpq(m,1/2-m));
   my(c = (-1)^m*b[1]/8 + S(m)*u[m+1]);
   print("  m=",m,"  ",valuation(c/u[m+1]-z2,2),"   ",valuation(v[m+1]/u[m+1]-z2,2),"   ",valuation((v[m+1]-c)/u[m+1],2)));
}
quit;
