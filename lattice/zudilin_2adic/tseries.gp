read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
PR=200; z2 = Lp(2,triv,2,PR);
tn(n) = (2^(n+1)-2)*bernfrac(n);
Tser(xx,NT) = sum(n=0,NT, tn(n)*(-1/xx)^(n+1));
NT=260;
print("T(1/2) v2 = ",valuation(Tser(1/2,NT),2));
print("T(1/2)/z2 = ", Tser(1/2,NT)/z2 + O(2^40));
print("T(1/2)+8z2 v2 = ", valuation(Tser(1/2,NT)+8*z2,2));
print("T(1/2)-8z2 v2 = ", valuation(Tser(1/2,NT)-8*z2,2));
print("shift: T(x)+T(x+1)-2/x^2 at x=-1/2: v2=", valuation(Tser(-1/2,NT)+Tser(1/2,NT)-2/(-1/2)^2,2));
print("shift: T(x)+T(x+1)+2/x^2 at x=-1/2: v2=", valuation(Tser(-1/2,NT)+Tser(1/2,NT)+2/(-1/2)^2,2));
bpqall(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
{for(m=0,5, my(b=bpqall(45,1/2-m), r=b[1][46]/b[2][46]);
  print("x=",1/2-m,"  v2(CFlimit - Tser)=",valuation(r-Tser(1/2-m,NT),2), "  v2(CFlimit + Tser)=",valuation(r+Tser(1/2-m,NT),2), "  v2(CF)=",valuation(r,2)));}
quit;
