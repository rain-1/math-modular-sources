default(realprecision,60);
s2=sqrt(2); s3=sqrt(3);
print("min polys of K pi^{3/2}:");
{foreach([[2,"alpha"],[sqrt((24+17*s2)/32),"gamma"],[sqrt((4+3*s2)/4),"eps"],[sqrt(9*(2+s3)/16),"zeta"],[sqrt(27/16),"s7"],[s2,"s10"],[3*s2,"s18"]],P,
  print("  ",P[2],": ",P[1],"  algdep4 = ",algdep(precision(P[1],45),4)));}
print("x_+ eps  = ",1/(12+8*s2),"  (3-2 sqrt2)/4 = ",(3-2*s2)/4);
print("x_+ zeta = ",1/(9+6*s3),"  (2 sqrt3-3)/9 = ",(2*s3-3)/9);
M=40;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
u6 = 'q*(et(2)*et(6)/(et(1)*et(3)))^6;
print("(eta2 eta6/(eta1 eta3))^6 = ", u6+O('q^8));
u12 = 'q*(et(2)*et(3)*et(12)/(et(1)*et(4)*et(6)))^4;
x12 = u12/(1+10*u12+9*u12^2);
print("alpha x(q) = ", x12+O('q^8));
print("alpha x(-q) + u6-quotient : ", subst(truncate(x12+O('q^9)),'q,-'q) + subst(truncate(u6+O('q^9)),'q,'q));
quit;
