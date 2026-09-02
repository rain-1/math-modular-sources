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
NN=Pol([2985984,-99035136,-80870400,165680640,306913536,216563328,82569024,18046944,2131344,95880,-3900,-398,1],'z);
DD=(2*'z+1)^4*(3*'z+1)^2*(4*'z+1)^2*(6*'z+1)^4;
lhs=R*subst(DD,'z,v)-subst(NN,'z,v);
print("residual R*Den(v) - Num(v) = ",lhs);
print("valuation of residual (should be >= 690): ",if(lhs==0,"identically 0 to series precision",valuation(lhs,'q)));
\\ kernel dimension at degree 12
fitk(GG,uu,d1,d2,NC)=
{
  my(cols=List(),Ms);
  for(j=0,d1, listput(cols, Vec(truncate(uu^j+O('q^(NC+1))),-(NC+1))));
  for(j=0,d2, listput(cols, Vec(truncate(-GG*uu^j+O('q^(NC+1))),-(NC+1))));
  Ms=Mat(vector(#cols,i,cols[i]~));
  matker(Ms);
}
print("kernel dim at (12,12) with 169 equations: ",#fitk(R,v,12,12,168));
print("kernel dim at (11,11): ",#fitk(R,v,11,11,164));
quit;
