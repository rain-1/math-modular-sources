default(parisize,"8G");
default(realprecision,120);
s2=sqrt(2); s3=sqrt(3);
W3=gamma(1/3)^6/Pi^4;   \\ weight-2 CM unit for Q(sqrt-3)
W4=gamma(1/4)^4/Pi^3;   \\ weight-2 CM unit for Q(i)
print("=== elementary rows: K = sqrt(N)/(2 pi^{3/2}) sqrt(lam1/(lam1-lam2))");
RS=[["alpha",12,16,4],["gamma",6,17+12*s2,17-12*s2],["eps",8,12+8*s2,12-8*s2],["zeta",9,9+6*s3,9-6*s3],["s7",7,27,-1],["s10",10,16,-4],["s18",18,16,12]];
CF=[2, (1+s2)^2/2^(9/4), sqrt(4+3*s2)/2, 3*(sqrt(2)+sqrt(6))/8, 3*sqrt(3)/4, s2, 3*s2];
{for(i=1,#RS, my(R=RS[i],N=R[2],l1=R[3],l2=R[4]);
  my(K=sqrt(N)/(2*Pi^(3/2))*sqrt(l1/(l1-l2)));
  print(R[1],":  K = ",K,"   K pi^{3/2} = ",K*Pi^(3/2),"   K^2 pi^3 = ",K^2*Pi^3);
  print("      closed form check: ", K*Pi^(3/2)-CF[i]));}
print();
print("=== level-12 zeta(5) [CM]");
Dv3 = (3*(45+26*s3)/2^17);      \\ (Dv/W3)^3
K5 = 9*(Dv3)^(2/3)*sqrt((7+4*s3)/(24*s3))*W3^2/sqrt(Pi);
print("  K = ", K5);
print("  K pi^{17/2}/Gamma(1/3)^12 = ", K5*Pi^(17/2)/gamma(1/3)^12);
print("  published K               = 0.68537184849053440595143004488054241747925748918627");
print();
print("=== level-16 zeta(5) [weight 4, Fricke-even, F(tau_c) != 0]");
Rc = (10021-7083*s2)/2;   \\ |R(x_c)|
K16 = Rc*(2+s2)*sqrt(4+s2)/8*W4/Pi^(3/2);
print("  |R(x_c)| = ", Rc);
print("  K = ", K16);
print("  K pi^{9/2}/Gamma(1/4)^4 = ", K16*Pi^(9/2)/gamma(1/4)^4);
print("  (numerical fold value was 2.0499712563401062935751009280207553319602675007961708038442)");
print();
print("=== level-12 zeta(7) [weight 6, Fricke-odd, F(tau_c) != 0]");
A1 = 2025*(849969-490730*s3);
A2 = (1152*s3-2016)/25;
K7 = (3*s3/2)*A1*(Dv3)^(2/3)*sqrt(1/(5*(-A2)))*W3^2/Pi^(3/2);
print("  A1 = F/(Dv)^3 = ", A1);
print("  A2 = D2t/(Dv)^2 = ", A2);
print("  K = ", K7);
print("  K pi^{19/2}/Gamma(1/3)^12 = ", K7*Pi^(19/2)/gamma(1/3)^12);
print("  (numerical fold value was 31.721520827964090772700383570287508787127452058803917818394)");
quit;
