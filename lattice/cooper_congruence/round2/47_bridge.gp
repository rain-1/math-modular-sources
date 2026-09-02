/* 47_bridge.gp -- consolidating checks.
   (a) the character formula   psi(p) = kappa^{-1} * A_{p-1}/p  (mod p),  kappa=-2,-3,-3
   (b) the bridge  kappa * a_{p-1} = A_{p-1}/p  (mod p)
   (c) deg_x F_{<p} mod p  vs  floor( mu(N)(p-1) / (6 deg x) )  -- the supersingular count
       (mu(7)=8, mu(10)=18, mu(18)=36;  deg x = 2,4,8)
   (d) the "b" form of the target restated with F_{<p}. */
default(parisize,10000000000);
read("40_core.gp");
NX = 1200;
KAP = [-2,-3,-3];
MU  = [8,18,36];
DGX = [2,4,8];
{ AAv=vector(3); EEv=vector(3);
  for(k=1,3, AAv[k]=AVEC(k,NX); EEv[k]=AJ(k,NX)); }
print("=== (a)/(b)  kappa*a_{p-1} = A_{p-1}/p = kappa*psi(p)  (mod p) ===");
{ for(k=1,3, my(A=AAv[k], a=EEv[k], f1=[], f2=[]);
   forprime(p=5,300, if(p>NX,break);
     my(ps=psival(k,p), t);
     if(A[p]%p!=0, f1=concat(f1,[p]); next);
     t = A[p]/p;
     if(Mod(t,p) != Mod(KAP[k]*ps,p), f1=concat(f1,[p]));
     if(Mod(KAP[k]*a[p],p) != Mod(t,p), f2=concat(f2,[p])));
   print("  ",NAM[k],":  A_{p-1}/p = kappa*psi(p) mod p, failures 5<=p<=300: ",if(#f1==0,"NONE",f1));
   print("        kappa*a_{p-1} = A_{p-1}/p mod p,  failures: ",if(#f2==0,"NONE",f2))); }
print();
print("=== (c) deg_x F_{<p} (mod p) vs the supersingular count floor(mu(N)(p-1)/(6 deg x)) ===");
{ for(k=1,3, my(A=AAv[k], row=[]);
   forprime(p=3,79, if(LEV[k]%p==0,next);
     my(Fp=sum(n=0,p-1, Mod(A[n+1],p)*'x^n), d=poldegree(Fp), pred=(MU[k]*(p-1))\(6*DGX[k]));
     row=concat(row,[[p,d,pred,if(d==pred,"=","DIFF")]]));
   print("  ",NAM[k]," [p, deg F_{<p}, predicted, match]: ",row)); }
print();
print("=== (d) restated target with F_{<p}:  first coefficients of H = P^((p-1)/2)/F_{<p} ===");
{ my(k=1, A=AAv[1], R=ROWS[1], B=R[3], C=R[4]);
  forprime(p=5,13,
    my(Fp=sum(n=0,p-1, Mod(A[n+1],p)*'x^n), Pw=Mod(1,p)*(1-2*B*'x+(B^2-4*C)*'x^2)^((p-1)/2), g);
    g = gcd(Pw,Fp);
    print("  s7 p=",p,":  P^((p-1)/2)=",lift(Pw/g),"   F_{<p}=",lift(Fp/g),"   (common factor deg ",poldegree(g),")")); }
quit;
