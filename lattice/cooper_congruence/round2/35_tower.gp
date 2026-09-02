\\ 35_tower.gp -- Task (3): numerical verification of the tower mechanism.
\\   (T1)  g_{m0}   |T_{p^2} = psi(p) p g_{m0} + p^3 g_{m0 p^2}
\\   (T2)  g_{m0p^{2r}}|T_{p^2} = g_{m0 p^{2r-2}} + p^3 g_{m0 p^{2r+2}}   (r>=1)
\\ Consequences to test on the lift side (weight 2k):
\\   (C1)  A(pn)   = lam A(n)     mod p^{2k-1}
\\   (C2)  A(p^2n) = lam A(pn)    mod p^{2(2k-1)}            <-- the extra cancellation
\\   (C3)  A(p^r n)= lam^r A(n)   mod p^{r+2}   (k=2) / p^{2r+3} (k=3)
\\   (C4)  psi(p)=1  ==>  v_p(a(p^{2j}m^2)) >= 3j + v_p(a(m^2))
default(parisize, 6000000000);
NN = 1300; MM = 1250;
{ E4s = 1 + 240*sum(n=1, NN-1, sigma(n,3)*'q^n) + O('q^NN); }
{ E6s = 1 - 504*sum(n=1, NN-1, sigma(n,5)*'q^n) + O('q^NN); }
DEL = (E4s^3 - E6s^2)/1728;
A4a = vector(MM,m,polcoeff(DEL/E4s^2,m));
A4b = vector(MM,m,polcoeff(E4s*DEL/E6s^2,m));
A6  = vector(MM,m,polcoeff(E6s*DEL/E4s^3,m));
{ Ac(V,n) = if(n<1 || n>#V, 0, V[n]); }
{ tst(nam, V, D0, kk, cden) = my(w2=2*kk-1);
  print("");
  print("=== ", nam, "  (weight ", 2*kk, ", lam = chi_",D0,"(p) p^{k-1} = psi(p) p^",kk-1,") ===");
  forprime(p=2,19,
    my(lam, b1=0,b2=0,b3=0,f1=0,f2=0, NT);
    lam = kronecker(D0,p)*p^(kk-1);
    NT = MM\p^2;
    for(n=1, MM\p, if((Ac(V,p*n) - lam*Ac(V,n)) % p^w2 != 0, b1++; if(f1==0,f1=n)));
    for(n=1, NT,   if((Ac(V,p^2*n) - lam*Ac(V,p*n)) % p^(2*w2) != 0, b2++; if(f2==0,f2=n)));
    for(r=1,4, for(n=1, MM\p^r,
      if((Ac(V,p^r*n) - lam^r*Ac(V,n)) % p^(r*(kk-1)+w2) != 0, b3++)));
    print("  p=",p," lam=",lam,
      "   (C1) A(pn)=lam A(n) mod p^",w2,": ", if(b1==0,"PASS",Str("FAIL ",b1," first ",f1)),
      "   (C2) A(p^2n)=lam A(pn) mod p^",2*w2,": ", if(b2==0,"PASS",Str("FAIL ",b2," first ",f2)),
      "   (C3) r<=4: ", if(b3==0,"PASS",Str("FAIL ",b3))));
}
tst("F4a = Delta/E4^2 (pole disc -3)", A4a, -3, 2, 64);
tst("F4b = E4 Delta/E6^2 (pole disc -4)", A4b, -4, 2, 108);
tst("F6 = E6 Delta/E4^3 (pole disc -3, k=3)", A6, 1, 3, 384);
print("");
print("=== (C4)  psi(p)=1  ==>  p^{3j} | a(p^{2j} m^2),  i.e.  v_p(beta(p^j m)) >= 2j ===");
{ bet(V) = my(M=#V, cp=vector(M,n,V[n]/n), b=vector(M));
  for(n=1,M, my(s=0); fordiv(n,d, s+=moebius(d)*cp[n/d]); b[n]=s); b; }
B4a = bet(A4a);
B4b = bet(A4b);
{ c4(nam, B, D0) = print("  --- ", nam, " ---");
  forprime(p=2,29, my(ps=kronecker(D0,p), bad=0, cnt=0, mn=99);
    for(n=p, #B, my(e=valuation(n,p)); if(e>0, cnt++;
      my(vv=if(B[n]==0,10^6,valuation(B[n],p)));
      if(vv < 2*e, bad++);
      if(vv-2*e<mn, mn=vv-2*e)));
    print("   p=",p," psi(p)=",ps," : v_p(beta(n))>=2 v_p(n) ", if(bad==0,Str("PASS (min excess ",mn,")"),Str("FAIL ",bad,"/",cnt))));
}
c4("F4a, psi=chi_{-3}", B4a, -3);
c4("F4b, psi=chi_{-4}", B4b, -4);
quit;
