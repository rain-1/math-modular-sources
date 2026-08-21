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
S(m) = sum(j=1,m,(-1)^j/(2*j-1)^2);
NN=40;
{ for(m=0,6,
   my(xx=1/2-m, b=bpqall(NN,xx), th, cauchy);
   th = b[1][NN+1]/b[2][NN+1];
   cauchy = valuation(th - b[1][NN]/b[2][NN], 2);
   print("x=",xx,"  Theta_approx cauchy v2=",cauchy);
   print("    v2(Th + 8*z2)=",valuation(th+8*z2,2), "   v2(Th - (-1)^m*(-8*z2-8*S(m)))=", valuation(th - (-1)^m*(-8*z2-8*S(m)),2), "  v2(Th - (-1)^m*(-8*z2+8*S(m)))=",valuation(th-(-1)^m*(-8*z2+8*S(m)),2));
  );
}
quit;
