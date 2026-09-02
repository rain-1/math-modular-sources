default(parisize,"20G");
default(realprecision,700);
cn=read("cn400.txt");
dn=read("dn.txt");
xi=(209/1728)*zeta(7);
NN=400;
e=vector(NN+1,i,dn[i]-xi*cn[i]);
f=vector(NN,i,e[i+1]+e[i]);
Lg=log(1+x+O(x^(NN+2)));
mkb(II)=
{
  my(bas=List(),lab=List());
  for(i=0,II, my(pf=(1+x+O(x^(NN+2)))^i); for(j=1,7, listput(bas,pf*Lg^j); listput(lab,[i,j])));
  [Vec(bas),Vec(lab)];
}
run(II,nodes)=
{
  my(bl=mkb(II), bas=bl[1], lab=bl[2], nb=#bas, M, rhs, so);
  if(#nodes!=nb, return("size mismatch"));
  M=matrix(nb,nb,r,s,polcoeff(bas[s],nodes[r]));
  rhs=vector(nb,r,f[nodes[r]])~;
  so=matsolve(M,rhs);
  for(k=1,nb, if(lab[k]==[0,7], return(so[k])));
}
geo(nb,lo,hi)=vector(nb,k,round(lo*(hi/lo)^((k-1)/(nb-1))));
{foreach([1,2,3,4],II, my(nb=7*(II+1));
 foreach([80,120,160,200],lo, my(nd=geo(nb,lo,400));
  if(#Set(nd)==nb, print("I=",II," nodes ",lo,"..400 : alpha = ",precision(run(II,nd),25)),
    print("I=",II," nodes ",lo,"..400 : duplicate nodes"))));}
print();
print("14161/5040 = ",precision(14161/5040.,25));
print("uniform-window checks, I=2:");
{foreach([200,240,280,320,359],n0, my(nd=vector(21,k,n0+(k-1)*2)); if(nd[21]<=NN, print("  n0=",n0," alpha=",precision(run(2,nd),25))));}
quit;
