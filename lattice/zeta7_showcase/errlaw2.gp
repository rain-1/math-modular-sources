default(parisize,"20G");
default(realprecision,700);
cn=read("cn400.txt");
dn=read("dn.txt");
xi=(209/1728)*zeta(7);
NN=400;
e=vector(NN+1,i,dn[i]-xi*cn[i]);
Lg=log(1+x+O(x^(NN+2)));
inv=1/(1+x+O(x^(NN+2)));
bas=List();
lab=List();
listput(bas,inv); listput(lab,[0,0]);
II=8;
{for(i=0,II, my(pref=inv*(1+x+O(x^(NN+2)))^i); for(j=1,7, listput(bas,pref*Lg^j); listput(lab,[i,j])));}
nb=#bas;
print("number of basis functions: ",nb);
B=vector(nb,k,vector(NN+1,m,polcoeff(bas[k],m-1)));
fit(n0)=
{
  my(M=matrix(nb,nb,r,s,B[s][n0+r]), rhs=vector(nb,r,e[n0+r])~);
  matsolve(M,rhs);
}
{foreach([NN+1-nb-40,NN+1-nb-20,NN+1-nb,NN+1-nb+1],n0,
  my(so=fit(n0));
  print("window n=",n0-1,"..",n0+nb-2,"  alpha=beta[0,7] = ",precision(so[nb],30)));}
print();
print("candidate 14161/5040 = ",precision(14161/5040.,30));
print("candidate 39350529/5040 = ",precision(39350529/5040.,20));
so=fit(NN+1-nb);
print("full solution vector (label -> coeff):");
{for(k=1,nb, print("  ",lab[k]," -> ",precision(so[k],25)));}
quit;
