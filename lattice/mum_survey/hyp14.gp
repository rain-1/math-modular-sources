/* The 14 hypergeometric CY operators: factorial-ratio data, the constant
   C = 1/t_c, slope primes, and the Gamma-class constant chi/kappa at every
   place (archimedean: zeta(3); p-adic: zeta_p(3) := L_p(3,omega^-2)).      */
read("gammap.gp");
read("../euler_criterion/lp.gp");
Lpom(p,m,s,PR)=
{ my(q=if(p==2,4,p),F=q,J,tot,br,inn,x,tt,c);
  J=ceil((PR+5)/valuation(F,p))+3; setbern(J); tot=O(p^PR);
  for(a=1,F, if(a%p==0,next);
    c=if(p==2, if(a%4==1,1,-1)^m, teichmuller(a+O(p^(PR+8)))^m);
    br=(a+O(p^(PR+8)))/tch(a,p,PR+8); x=(F/a)+O(p^(PR+8));
    inn=O(p^(PR+8)); tt=1+O(p^(PR+8));
    for(j=0,J, inn+=binomial(1-s,j)*tt*bern(j); tt*=x);
    tot+=c*br^(1-s)*inn);
  tot/(F*(s-1)); }

/* primitive set at N contributes  {(d, mu(N/d))} */
prim(N) = my(L=List()); fordiv(N,d, my(m=moebius(N/d)); if(m, listput(L,[d,m]))); Vec(L);
merge(lists) = { my(H=Map()); for(i=1,#lists, my(L=lists[i]);
    for(j=1,#L, my(u=L[j][1], e=L[j][2], cur=0);
      if(mapisdefined(H,u), cur=mapget(H,u)); mapput(H,u,cur+e)));
  my(ks=Set(Vec(Mat(H)[,1])), R=List());
  for(i=1,#ks, my(u=ks[i], e=mapget(H,u)); if(e, listput(R,[u,e])));
  vecsort(Vec(R),1); }

/* the 14 weight systems, as lists of primitive sets (with multiplicity) */
{FAM = [
 ["1/5,2/5,3/5,4/5",   [5]],
 ["1/8,3/8,5/8,7/8",   [8]],
 ["1/10,3/10,7/10,9/10",[10]],
 ["1/12,5/12,7/12,11/12",[12]],
 ["1/2,1/2,1/2,1/2",   [2,2,2,2]],
 ["1/3,1/3,2/3,2/3",   [3,3]],
 ["1/4,1/4,3/4,3/4",   [4,4]],
 ["1/6,1/6,5/6,5/6",   [6,6]],
 ["1/3,1/2,1/2,2/3",   [3,2,2]],
 ["1/4,1/2,1/2,3/4",   [4,2,2]],
 ["1/6,1/2,1/2,5/6",   [6,2,2]],
 ["1/4,1/3,2/3,3/4",   [4,3]],
 ["1/6,1/3,2/3,5/6",   [6,3]],
  ["1/6,1/4,3/4,5/6",   [6,4]]
];}
{
print("weights | factorial exponents (u^e) | C=1/t_c | slope primes | sum e*u^3 | chi/kappa");
for(i=1,#FAM,
  my(f=FAM[i], Ns=f[2], LL=List());
  for(j=1,#Ns, listput(LL, prim(Ns[j])));
  listput(LL, [[1,-4]]);
  my(E=merge(Vec(LL)), C=1, s3=0, str="");
  for(j=1,#E, my(u=E[j][1], e=E[j][2]);
    C *= u^(e*u); s3 += e*u^3;
    str = Str(str, if(j>1,"*",""), "(",u,"n)!^",e));
  print(f[1]," | ",str," | ",C," | ",factor(C)[,1]~," | ",s3," | ",-s3/3);
);
}
print();
print("=== p-adic check: [rho^3] log(Gamma_p-class) = (chi/kappa) * zeta_p(3) ? ===");
{ my(PR=200,K=40);
  for(pi=1,4, my(p=[2,3,5,7][pi], c=gpTaylor(p,K,PR), z3=Lpom(p,-2,3,60), c3=c[3]);
    print("p=",p,"   v_p(c3 + zeta_p(3)/3) = ", valuation(c3 + z3/3, p),
          "   [c3 = -zeta_p(3)/3 <=> sum e u^3 * c3 = (chi/kappa) zeta_p(3)]");
  );
}
quit
