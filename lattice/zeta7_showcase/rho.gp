default(parisize,"32G");
NQ=700;
et(k)=prod(m=1,NQ\k+1,1-'q^(k*m)+O('q^(NQ+2)));
D(s)='q*deriv(s,'q);
h=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/'q;
v=1/h;
x=h/((h+3)*(h+4));
E8(k)=1+480*sum(m=1,NQ\k+1,sigma(m,7)*'q^(k*m))+O('q^(NQ+2));
Phi=(E8(1)-572*E8(2)+11583*E8(3)-36608*E8(4)+46332*E8(6)-20736*E8(12))/480;
U=Phi/D(x);
R=v^3*U/D(v)^3;
print("R = v^3 U/(Dv)^3 = ",R+O('q^8));
fit(GG,uu,d1,d2,NC)=
{
  my(cols=List(),Ms);
  for(j=0,d1, listput(cols, Vec(truncate(uu^j+O('q^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-GG*uu^j+O('q^(NC+1))),-(NC+1))));
  Ms=Mat(vector(#cols,i,cols[i]~));
  matker(Ms);
}
{for(d=1,16, my(NC=4*d+120, k=fit(R,v,d,d,NC));
  if(#k>0, print("R rational in v of degree ",d,":");
    print("  P = ",vector(d+1,j,k[j,1]));
    print("  Q = ",vector(d+1,j,k[d+1+j,1]));
    break));}
quit;
