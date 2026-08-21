\\ homogeneous fit: find (N,D) with deg N<=dn, deg D<=dd, N(x) - g*D(x) = 0
ratfit(g,x,dn,dd,nc)={
 my(vars=dn+dd+2, M, xp=vector(dn+dd+1));
 xp[1]=1+O(q^PREC); for(i=1,dn+dd, xp[i+1]=xp[i]*x);
 M=matrix(nc,vars);
 for(m=0,nc-1,
   for(i=0,dn, M[m+1,i+1]=polcoeff(xp[i+1],m));
   for(j=0,dd, M[m+1,dn+2+j]=-polcoeff(g*xp[j+1],m)));
 my(K=matker(M));
 if(#K==0, return(0));
 my(s=K[,1], N=sum(i=0,dn,s[i+1]*x^i), D=sum(j=0,dd,s[dn+2+j]*x^j));
 if(#K>1, return([-1,#K]));
 my(err=N-g*D); if(err!=0 && valuation(err,q)<nc-2, return(0));
 [N,D,#K] };
scan(g,x,nm,mx,nc)={my(found=0); for(dd=0,mx, for(dn=0,mx, my(r=ratfit(g,x,dn,dd,nc)); if(type(r)=="t_VEC" && r[1]!=-1, print(nm,": dn=",dn," dd=",dd,"  N=",r[1],"  D=",r[2]); found=1; break(2)))); found}
