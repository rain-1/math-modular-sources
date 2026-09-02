default(parisize,"12G");
M = 260;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
D(s)='q*deriv(s,'q);
E4(k)=1+240*sum(m=1,M\k+1,sigma(m,3)*'q^(k*m))+O('q^(M+2));
E6(k)=1-504*sum(m=1,M\k+1,sigma(m,5)*'q^(k*m))+O('q^(M+2));
E8(k)=1+480*sum(m=1,M\k+1,sigma(m,7)*'q^(k*m))+O('q^(M+2));
rf(GG,U,d1,d2,NC)={
  my(cols=List(), Ms);
  for(j=0,d1, listput(cols, Vec(truncate(U^j+O('q^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-GG*U^j+O('q^(NC+1))),-(NC+1))));
  Ms = Mat(vector(#cols,i,cols[i]~));
  matker(Ms);
}
show(nm,RR,U,dmax)={
  print("--- ",nm);
  print("    R = ", RR+O('q^7));
  for(d=0,dmax, my(k=rf(RR,U,d,d,6*d+60)); if(#k>0, print("    fit degree ",d,", kernel dim ",#k,": ", k~); break));
}
\\ level 16 zeta(5)
X16 = 'q*et(2)*et(16)^2/(et(1)^2*et(8));
T16 = X16/(8*X16^2+2*X16+1);
PH16 = -(E6(1)-85*E6(2)+1428*E6(4)-5440*E6(8)+4096*E6(16))/504;
F16 = PH16/D(T16);
show("level16 zeta5:  R = x^2 F/(Dx)^2", X16^2*F16/D(X16)^2, X16, 8);
\\ level 12 zeta(7)
X12 = 'q*et(4)^2*et(12)^2/(et(1)^2*et(3)^2);
T12 = X12/(16*X12^2+2*X12+1);
PH12 = (E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
F12 = PH12/D(T12);
show("level12 zeta7:  R = x^3 F/(Dx)^3", X12^3*F12/D(X12)^3, X12, 8);
\\ level 12 zeta(5)
h = et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/'q;
v = 1/h;
Xz5 = v/(1+7*v+12*v^2);
E = (-9*E4(3)+16*E4(4))/7;
G = 7*E/(143*Xz5^3+189*Xz5^2+21*Xz5+7);
show("level12 zeta5:  R = v^2 G/(Dv)^2", v^2*G/D(v)^2, v, 8);
quit;
