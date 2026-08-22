\\ exact rho via Q(sqrt(-3)):  w = sqrt(-3) = Mod(y, y^2+3)
{
  my(w = Mod(y, y^2+3));
  \\ T1
  my(a=117,b=78,c=21,d=3969,e=-2646,f=441);
  my(tc = (13 - 3*w)/882);
  print("T1 exact: R(tc) = ", lift(1 - a*tc + d*tc^2));
  print("T1 exact: tc*R'(tc) = ", lift(tc*(-a+2*d*tc)));
  print("T1 exact: T(tc) = ", lift(1 - b*tc + (d+e)*tc^2));
  print("T1 exact: rho = ", lift(-(1 - b*tc + (d+e)*tc^2)/(tc*(-a+2*d*tc))));
  print("T1 exact: |tc|^2 = ", lift(tc*(13+3*w)/882));
  \\ T2
  a=72;b=36;c=6;d=1728;e=-1728;f=324;
  tc = (3 - w)/144;
  print("T2 exact: R(tc) = ", lift(1 - a*tc + d*tc^2));
  print("T2 exact: tc*R'(tc) = ", lift(tc*(-a+2*d*tc)));
  print("T2 exact: T(tc) = ", lift(1 - b*tc + (d+e)*tc^2));
  print("T2 exact: rho = ", lift(-(1 - b*tc + (d+e)*tc^2)/(tc*(-a+2*d*tc))));
  print("T2 exact: |tc|^2 = ", lift(tc*(3+w)/144));
}
s3 = sqrt(3);
tc1 = (13 - 3*I*s3)/882;
xi1 = dofold("T1 Herfurtner #30 (Gamma_0(7))", 117,78,21, 3969,-2646,441, tc1, 2500, 1200, 2400);
tc2 = (3 - I*s3)/144;
xi2 = dofold("T2 Herfurtner #45 (level 3)", 72,36,6, 1728,-1728,324, tc2, 2500, 1100, 2400);
write("xi_targets.txt", "T1 " , xi1);
write("xi_targets.txt", "T2 " , xi2);
