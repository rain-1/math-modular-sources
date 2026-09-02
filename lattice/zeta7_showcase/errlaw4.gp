default(parisize,"40G");
default(realprecision,1150);
cn=read("cn.txt");
dn=read("dn900.txt");
xi=(209/1728)*zeta(7);
NN=900;
e=vector(NN+1,i,dn[i]-xi*cn[i]);
f=vector(NN,i,e[i+1]+e[i]);
Lg=log(1+x+O(x^(NN+2)));
mkb(II)=
{
  my(bas=List(),lab=List());
  for(i=0,II, my(pf=(1+x+O(x^(NN+2)))^i); for(j=1,7, listput(bas,pf*Lg^j); listput(lab,[i,j])));
  [Vec(bas),Vec(lab)];
}
BL=[mkb(1),mkb(2),mkb(3),mkb(4),mkb(5)];
run(II,nodes)=
{
  my(bl=BL[II], bas=bl[1], lab=bl[2], nb=#bas, M, rhs, so);
  M=matrix(nb,nb,r,s,polcoeff(bas[s],nodes[r]));
  rhs=vector(nb,r,f[nodes[r]])~;
  so=matsolve(M,rhs);
  for(k=1,nb, if(lab[k]==[0,7], return(so[k])));
}
geo(nb,lo,hi)=vector(nb,k,round(lo*(hi/lo)^((k-1)/(nb-1))));
print("target 2023/720 = ",precision(2023/720.,30));
{foreach([1,2,3,4,5],II, my(nb=7*(II+1));
 foreach([100,200,300,450],lo, my(nd=geo(nb,lo,900));
  if(#Set(nd)==nb, print("I=",II-1," nodes ",lo,"..900 : alpha = ",precision(run(II,nd),30)))));}
print();
print("raw ratios (-1)^(n+1) f_n / [x^n]log^7(1+x):");
L7=Lg^7;
{foreach([100,200,300,400,500,600,700,800,900],n, print("  n=",n,"  ",precision(f[n]/polcoeff(L7,n),20)));}
quit;
