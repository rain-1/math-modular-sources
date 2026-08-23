/* run_h5.gp -- Task 4 part 2: do the *targets* converge?
   (a) the 2-adic congruence target  eta_n = -U_n/Y_n  (unit part), and
       -V_n/X_n; their agreement is the 2-adic bridge;
   (b) the archimedean target  rho_n = -(V_nG-U_n)/(X_nG-Y_n);
   (c) is xi2 = zeta_2(2) algebraic 2-adically?  (Ridout/p-adic Roth need it) */
\p 3000
default(parisize, 6000000000);
GG = Catalan;
DIR = "/home/ubuntu/code/math-modular-sources/lattice/p2_holonomic/data/";
RW  = rdrows(concat(DIR,"rows_n200.txt"));
XD = xi2at(400); XI = XD[1]; XR = XD[3];
print("--- is xi2 = zeta_2(2) 2-adically algebraic?  algdep at 2-adic prec 600 ---");
{
my(z = XR + O(2^600));
for(d=1,8, my(pp=algdep(z,d));
  print("  deg ",d,": ",pp,"   log10 height = ",
        if(poldegree(pp)>0, round(log(vecmax(abs(Vec(pp))))/log(10)), 0)));
}
print("");
print("--- lindep([1,xi2]) at 2-adic precision 600 ---");
print(lindep([1+O(2^600), XR+O(2^600)]));
print("");
print("n,v2etaU,v2etaV,v2diff,unitU40,unitV40,logrhoinf,v2XiForm,v2ViForm");
{
for(nn=4,200,
 my(rw=mapget(RW,nn), XX=rw[1], YY=rw[2], VV=rw[3], UU=rw[4],
    eU=-UU/YY, eV=-VV/XX, vU=valuation(eU,2), vV=valuation(eV,2),
    uU=lift((eU/2^vU)+O(2^40)), uV=lift((eV/2^vV)+O(2^40)),
    dif=valuation(eU-eV,2),
    rin=log((VV*GG-UU)/abs(XX*GG-YY)));
 printf("%d,%d,%d,%d,%d,%d,%.6f,%d,%d\n", nn, vU, vV, dif, uU, uV, rin,
    valuation(XX*XR-YY,2), valuation(VV*XR-UU,2)));
}
\q
