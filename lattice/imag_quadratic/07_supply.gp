\\ 07_supply.gp -- the SUPPLY side.  For every level N <= 60 and every nebentypus eps
\\ of order 3, 4 or 6 (the orders whose character values generate an imaginary
\\ quadratic field: Q(zeta_3) for 3 and 6, Q(i) for 4) list the Eisenstein directions
\\ E_k^{psi1,psi2}(d tau) in the eps-component and the dimension of the fold-regular
\\ part (a_0 = 0 at infinity AND at the cusp 0), together with the field of definition.
\\ Rules used (proved in the report, sec.4; PARI-verified in 08_verify.gp):
\\   a_0(infinity) != 0  <=>  psi1 = 1        a_0(cusp 0) != 0  <=>  psi2 = 1
\\   hence every "mixed" direction (psi1 != 1 != psi2) is fold-regular, and each of the
\\   two pure blocks loses exactly one dimension.
\\ weight k = 3 for eps odd (period L(2,psi)), k = 4 for eps even (period L(3,psi)).
{
print("N    k  ord(eps) cond  #outer #inner #mixed   dim(foldreg)   field  minimal-K direction?");
for(N=3,60,
  my(G=znstar(N,1));
  for(me=1,N-1, if(gcd(me,N)>1,next);
    my(chi=znconreylog(G,me), ord=charorder(G,chi));
    if(ord!=3 && ord!=4 && ord!=6, next);
    my(par=chareval(G,chi,-1), k=if(par==0,4,3));
    my(f=znconreyconductor(G,chi));
    if(type(f)=="t_VEC", f=f[1]);
    my(nout=numdiv(N/f), nin=nout, nmix=0, fieldmix=[], dmix=0);
    for(m1=1,N-1, if(gcd(m1,N)>1,next);
      my(c1=znconreylog(G,m1), m2=lift(Mod(me,N)/Mod(m1,N)), c2=znconreylog(G,m2));
      my(f1=znconreyconductor(G,c1), f2=znconreyconductor(G,c2));
      if(type(f1)=="t_VEC", f1=f1[1]); if(type(f2)=="t_VEC", f2=f2[1]);
      if(f1==1 || f2==1, next);
      if(N%(f1*f2)!=0, next);
      my(o=lcm(charorder(G,c1),charorder(G,c2)));
      nmix += numdiv(N/(f1*f2));
      if(o==3||o==4||o==6, dmix += numdiv(N/(f1*f2)); fieldmix=concat(fieldmix,[o]));
    );
    my(dfr = nmix + max(0,nout-1) + max(0,nin-1));
    my(dK  = dmix + if(nout>1,1,0) + if(nin>1,1,0));
    if(dfr>0,
      print(N,"   ",k,"    ",ord,"      ",f,"      ",nout,"      ",nin,"      ",nmix,
            "        ",dfr,"        ",if(ord==4,"Q(i)","Q(zeta_3)"),
            "   dim_K = ",dK));
  );
);
}
quit;
