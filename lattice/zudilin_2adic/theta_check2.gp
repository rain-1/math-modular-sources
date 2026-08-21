read("/home/ubuntu/code/math-modular-sources/lattice/euler_criterion/lp.gp");
bpqall(n,xx) =
{ my(P=vector(n+2), Q=vector(n+2));
  Q[1]=1; Q[2]=xx^2-xx+1; P[1]=0; P[2]=1;
  for(m=1,n,
    Q[m+2]=((2*m*(m+1)+1-xx+xx^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
    P[m+2]=((2*m*(m+1)+1-xx+xx^2)*P[m+1]-m^2*P[m])/(m+1)^2);
  [P,Q];
}
PR=300; z2 = Lp(2,triv,2,PR);
print("z2 = ",z2);
b=bpqall(40,1/2);
{for(n=1,12, print("n=",n," v2(q)=",valuation(b[2][n+1],2)," v2(p)=",valuation(b[1][n+1],2)," ratio v2 cauchy=",if(n>1,valuation(b[1][n+1]/b[2][n+1]-b[1][n]/b[2][n],2))));}
th=b[1][41]/b[2][41];
print("Theta(1/2) approx v2 = ",valuation(th,2));
print("th/z2 = ", th/z2 + O(2^40));
print("th = ", th+O(2^30));
print("-8*z2 = ", -8*z2+O(2^30));
quit;
