/* 46_hexplicit.gp -- H made completely explicit, and the Dwork tower.
   The Lucas property A_{pj+r} = A_j A_r (mod p) (verified 0 failures, all rows, 44_lucas.log;
   Malik-Straub) says  F(x) = F_{<p}(x) F(x^p) (mod p),  F_{<p} := sum_{n<p} A_n x^n.
   Hence  F^{p-1} = F(x^p)/F(x) = 1/F_{<p}(x)  and
        H := (sqrt(P) F)^{p-1} = P^((p-1)/2) / F_{<p}(x)   (mod p)   -- a RATIONAL FUNCTION
   with numerator and denominator of degree < p.  The target congruence then reads
        [x^{pn}] ( l(x) * P(x)^((p-1)/2) / F_{<p}(x) ) = psi(p) L_n  (mod p),   n>=1.
   Also test the Dwork tower  a_{pn-1} = psi(p) a_{n-1}  (mod p^{1+v_p(n)}).            */
default(parisize,10000000000);
read("40_core.gp");
NX = 1500;
{ FFv=vector(3); LLv=vector(3); PPv=vector(3); AAv=vector(3); EEv=vector(3);
  for(k=1,3, AAv[k]=AVEC(k,NX); FFv[k]=FSER(k,NX); LLv[k]=LSER(k,NX); PPv[k]=PSER(k,NX);
             EEv[k]=LLv[k]/('x*sqrt(PPv[k])*FFv[k])); }
print("=== 1.  F^(p-1) = 1/F_{<p}  and  H = P^((p-1)/2)/F_{<p}  (mod p) ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], A=AAv[k], F=FFv[k], row=[]);
   forprime(p=3,23, if(LEV[k]%p==0,next);
     my(Fp = sum(n=0,p-1, Mod(A[n+1],p)*'x^n), T, ok=1, NT=400);
     T = (Mod(1,p)*(F+O('x^NT)))^(p-1) * Fp - 1;
     for(i=0,NT-2*p-4, if(polcoeff(T,i)!=0, ok=0; break));
     row=concat(row,[[p, poldegree(Fp), if(ok,"OK","FAIL")]]));
   print("  ",NAM[k]," [p, deg F_{<p} mod p, F^(p-1)*F_{<p}=1 ?]: ",row)); }
print();
print("=== 2.  the fully explicit target:  [x^{pn}]( l * P^((p-1)/2) / F_{<p} ) = psi(p) L_n ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], A=AAv[k], l=LLv[k], row=[]);
   forprime(p=3,53, if(LEV[k]%p==0,next);
     my(Fp = sum(n=0,p-1, Mod(A[n+1],p)*'x^n), Pw = Mod(1,p)*(1-2*B*'x+(B^2-4*C)*'x^2)^((p-1)/2),
        NT=1200, S, ps=psival(k,p), nf=0, nt=0);
     S = (Mod(1,p)*(l+O('x^NT))) * Pw / Fp;
     for(n=1,(NT-2*p-4)\p, nt++; if(polcoeff(S,p*n) != ps*Mod(polcoeff(l,n),p), nf++));
     row=concat(row,[[p,nt,nf]]));
   print("  ",NAM[k]," [p,#tests,#fail]: ",row)); }
print();
print("=== 3.  the Dwork tower:  a_{pn-1} = psi(p) a_{n-1}  (mod p^{1+v_p(n)}) ===");
{ for(k=1,3, my(e=EEv[k], ps, row=[]);
   forprime(p=2,31, my(nf=0, nt=0, worst=-1);
     ps = psival(k,p);
     for(n=1,(NX-2)\p, my(v=valuation(n,p), d=polcoeff(e,p*n-1)-ps*polcoeff(e,n-1));
       nt++;
       if(d!=0 && valuation(d,p) < 1+v, nf++; if(v>worst,worst=v)));
     row=concat(row,[[p,nt,nf]]));
   print("  ",NAM[k]," [p,#tests,#fail]: ",row)); }
print();
print("=== 4.  C(eta_1) = 0 mod p, i.e. a_{p^2 n -1} = psi(p) a_{pn-1} (mod p^2) ===");
{ for(k=1,3, my(e=EEv[k], row=[]);
   forprime(p=2,31, my(ps=psival(k,p), nf=0, nt=0);
     for(n=1,(NX-2)\(p^2), nt++;
       if(Mod(polcoeff(e,p^2*n-1), p^2) != Mod(ps*polcoeff(e,p*n-1),p^2), nf++));
     row=concat(row,[[p,nt,nf]]));
   print("  ",NAM[k]," [p,#tests,#fail]: ",row)); }
quit;
