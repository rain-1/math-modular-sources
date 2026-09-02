/* 13_haupt.gp -- numerical values of the Hauptmodul x at the four cusps of Gamma_1(5),
   confirming  x(oo)=0, x(0)=t1=phi^-5, x(1/2)=t2=-phi^5, x(2/5)=oo.        */
default(realprecision, 40);
leg5(n) = kronecker(n,5);
xval(tau, M) = my(qq = exp(2*Pi*I*tau)); qq*prod(n=1, M, (1-qq^n)^(5*leg5(n)));
s5 = sqrt(5); ph5 = (11+5*s5)/2; phm5 = (5*s5-11)/2;
print("t1 = phi^-5 = ", phm5, "   t2 = -phi^5 = ", -ph5);
print("");
print("towards cusp 0   (tau = i*y):");
for(k=1,5, my(y=[0.10,0.05,0.03,0.02,0.015][k]); print("   y=",y,"  x = ", xval(I*y, 4000)));
print("towards cusp 1/2 (tau = 1/2 + i*y):");
for(k=1,5, my(y=[0.10,0.05,0.03,0.02,0.015][k]); print("   y=",y,"  x = ", xval(1/2+I*y, 4000)));
print("towards cusp 2/5 (tau = 2/5 + i*y):");
for(k=1,4, my(y=[0.10,0.05,0.03,0.02][k]); print("   y=",y,"  x = ", xval(2/5+I*y, 4000)));
print("towards cusp oo  (tau = i*y):");
for(k=1,3, my(y=[1.0,2.0,3.0][k]); print("   y=",y,"  x = ", xval(I*y, 300)));
quit;
