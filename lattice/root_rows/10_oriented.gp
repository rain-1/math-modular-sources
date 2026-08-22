\\ Root rows in the *oriented* (elliptic) coordinate, where Sym^w really acts.
read("lattice/root_rows/lib.gp");
default(seriesprecision,400); default(realprecision,40);
M=400; q='q;
ee(d)={my(r=1+O(q^(M+1))); for(n=1,(M\d)+1, r*=(1-q^(d*n))); r};
E1=ee(1);E2=ee(2);E3=ee(3);E4=ee(4);E8=ee(8);
zag(a,b,c,N)={my(A=vector(N+1));A[1]=1;A[2]=b;for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2);A};

print("=========== beta(4) level 24 in the ORIENTED Catalan coordinate ===========");
{
 my(t = q*E1^4*E4^2*E8^4/E2^10, F = E2^10/(E1^4*E4^4));
 my(Qt=serreverse(t), Ft=subst(F,q,Qt));
 my(NB=300, g=vector(NB+1,j,polcoeff(Ft,j-1)));
 print("g (= cube root, in t) = ",vector(9,j,g[j]));
 print("integral: ",sum(i=1,NB+1,denominator(g[i])!=1)==0,"   lambda = 1");
 print("equals Zagier row E (12,4,32): ", g==vector(NB+1,j,zag(12,4,32,NB+1)[j]));
 my(mr=minrec(g,3,4)); print("root-row recurrence: order=",mr[1]," deg=",mr[2],"  ",mr[3]);
 print("root char roots: ",charroots(mr[3]));
 my(A3=vector(NB+1,j,polcoeff(Ft^3,j-1)));
 print("parent Sym^3 row in t: ",vector(9,j,A3[j]));
 print("parent integral: ",sum(i=1,NB+1,denominator(A3[i])!=1)==0);
 my(m3=minrec(A3,4,6)); print("parent recurrence: order=",m3[1]," deg=",m3[2]);
 print("parent char roots: ",charroots(m3[3]));
}
print("");
print("=========== L(4,chi_-3) level 24 in the ORIENTED level-3 coordinate ===========");
{
 my(t0 = q*(E3/E1)^12, b = E1^3/E3);
 my(Qt0=serreverse(t0), bt=subst(b,q,Qt0));
 my(NB=300, g=vector(NB+1,j,polcoeff(bt,j-1)));
 print("g (= cube root, in t0) = ",vector(9,j,g[j]));
 print("integral: ",sum(i=1,NB+1,denominator(g[i])!=1)==0,"   lambda = 1");
 my(mr=minrec(g,3,4)); print("root-row recurrence: order=",mr[1]," deg=",mr[2],"  ",mr[3]);
 print("root char roots: ",charroots(mr[3]));
 print("=> hypergeometric term: g_n = (-27)^n ((1/3)_n/n!)^2 ; check: ",
   vector(6,j,my(n=j-1); (-27)^n*(prod(i=0,n-1,(1/3+i))/n!)^2)==vector(6,j,g[j]));
 my(A3=vector(NB+1,j,polcoeff(bt^3,j-1)));
 my(m3=minrec(A3,4,6)); print("parent Sym^3 recurrence: order=",m3[1]," deg=",m3[2]);
 print("parent char roots: ",charroots(m3[3]));
}
\q
