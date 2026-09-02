default(parisize,"12G");
M = 400;
et(k)=prod(m=1,M\k+1,1-'q^(k*m)+O('q^(M+2)));
D(s)='q*deriv(s,'q);
E6(k)=1-504*sum(m=1,M\k+1,sigma(m,5)*'q^(k*m))+O('q^(M+2));
E8(k)=1+480*sum(m=1,M\k+1,sigma(m,7)*'q^(k*m))+O('q^(M+2));
rf(GG,U,d1,d2,NC)={
  my(cols=List(), Ms);
  for(j=0,d1, listput(cols, Vec(truncate(U^j+O('q^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-GG*U^j+O('q^(NC+1))),-(NC+1))));
  Ms = Mat(vector(#cols,i,cols[i]~));
  matker(Ms);
}
X16 = 'q*et(2)*et(16)^2/(et(1)^2*et(8));
T16 = X16/(8*X16^2+2*X16+1);
PH16 = -(E6(1)-85*E6(2)+1428*E6(4)-5440*E6(8)+4096*E6(16))/504;
F16 = PH16/D(T16);
R16 = X16^2*F16/D(X16)^2;
print("=== level16 zeta5");
{for(d=0,14, my(k=rf(R16,X16,d,d,4*d+120)); if(#k>0, print("  degree ",d," kernel dim ",#k,": ", k~); break));}
X12 = 'q*et(4)^2*et(12)^2/(et(1)^2*et(3)^2);
T12 = X12/(16*X12^2+2*X12+1);
PH12 = (E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
F12 = PH12/D(T12);
R12 = X12^3*F12/D(X12)^3;
print("=== level12 zeta7");
{for(d=0,14, my(k=rf(R12,X12,d,d,4*d+120)); if(#k>0, print("  degree ",d," kernel dim ",#k,": ", k~); break));}
quit;
