\\ 09_census.gp -- exact PARI census (no formula): for every N <= 30 and every
\\ nebentypus eps of order 3, 4 or 6 (character values in Q(zeta_3) or Q(i)) the
\\ dimension of the fold-regular part of the eps-component of M_k^Eis(Gamma_1(N)),
\\ k = 3 for eps odd (period L(2,psi)), k = 4 for eps even (period L(3,psi)).
\\ Fold-regular = a_0 vanishes at the cusp infinity AND at the cusp 0.
default(realprecision,60);
cycord(m) = {my(n,v=variable(m)); for(n=1,300, if(polcyclo(n,v)==m, return(n))); 0;}
toC(z) = {if(type(z)!="t_POLMOD", return(z*1.0)); my(m=z.mod, n=cycord(m)); subst(lift(z), variable(m), exp(2*Pi*I/n));}
print("N    k  ord(eps) cond(eps)  dim(eps-comp)  rank  dim(fold-regular)   K");
{
for(N=3,30,
  my(G=znstar(N,1));
  for(me=2,N-1, if(gcd(me,N)>1,next);
    my(chi=znconreylog(G,me), ord=charorder(G,chi));
    if(ord!=3 && ord!=4 && ord!=6, next);
    my(k=if(chareval(G,chi,-1)==0,4,3));
    my(f=znconreyconductor(G,chi));
    if(type(f)=="t_VEC", f=f[1]);
    my(mf=mfinit([N,k,[G,chi]],4), B=mfbasis(mf), d=#B, M, r);
    if(d==0, next);
    M=matrix(2,d);
    for(j=1,d, M[1,j]=toC(mfcoefs(B[j],0)[1]));
    for(j=1,d, my(p=0); M[2,j]=toC(mfslashexpansion(mf,B[j],[0,-1;1,0],0,1,&p)[1]));
    r=matrank(M);
    print(N,"    ",k,"     ",ord,"        ",f,"         ",d,"           ",r,"        ",d-r,"          ",if(ord==4,"Q(i)","Q(zeta_3)"));
  );
);
}
quit;
