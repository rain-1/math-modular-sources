default(realprecision,80);
K=0.68537184849053440595143004488054241747925748918627;
s3=sqrt(3); xp=7-4*s3; hc=2*s3;
Q=143*xp^3+189*xp^2+21*xp+7;
xpp=-2*hc/((hc+3)*(hc+4))^2;          \\ x''(h) at h^2=12
u=243*(23-15*s3)/14336;               \\ DE(tau*)/W^6
a=-(81*(45+26*s3)/2048)^(1/3);        \\ Dh(tau*)/W^2
W2=gamma(1/3)^6/Pi^4;
D2x=xpp*(a*W2)^2; DE=u*W2^3;
K2=49*xp*DE^2/(2*Pi*Q^2*(-D2x));
print("K (formula) = ",sqrt(K2));
print("K (numeric) = ",K);
print("ratio-1 = ",sqrt(K2)/K-1);
t=7*abs(u)*sqrt(xp/(-xpp))/(sqrt(2)*Q*abs(a));
print("t = K pi^{17/2}/Gamma(1/3)^12 = ",t,"   check: ",K*Pi^(17/2)/gamma(1/3)^12);
print("t minimal polynomial (algdep deg 12): ",algdep(t,12));
print("t^2 algdep deg 6: ",algdep(t^2,6));
print("t^6 algdep deg 2: ",algdep(t^6,2));
print("Q(x+) = ",Q,"  Q algdep: ",algdep(Q,2),"   xpp = ",xpp," algdep: ",algdep(xpp,2), "   xp/(-xpp) = ",xp/(-xpp)," algdep: ",algdep(xp/(-xpp),2));
quit;
