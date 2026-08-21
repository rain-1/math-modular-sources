\\ find linear relations among a list of q-series (as vectors of coefficients)
relker(L,nc)={ my(M=matrix(nc,#L)); for(j=1,#L, for(m=0,nc-1, M[m+1,j]=polcoeff(L[j],m))); matker(M) }
\\ monomials u^i v^j, i<=di, j<=dj
monos(u,v,di,dj)={my(L=List(),T=List()); for(i=0,di, for(j=0,dj, listput(L,u^i*v^j); listput(T,[i,j]))); [Vec(L),Vec(T)]}
