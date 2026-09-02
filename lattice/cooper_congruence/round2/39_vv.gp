\\ 39_vv.gp -- the vector-valued reading: two cheap decisive computations.
\\ (A) obstacle (i): is the relevant holomorphic space zero?
\\     S^+_{5/2}(Gamma_0(28)) = Shimura preimage of S_4(Gamma_0(7)); the s7 object is
\\     W_7-ANTIinvariant (Phi|W_7 = -Phi), so what matters is the eigenvalue of the
\\     weight-4 newform on Gamma_0(7) under the Fricke/Atkin-Lehner involution W_7.
\\ (B) the T_{p^2} test on the twisted-trace coefficients c(d) of row s7.
default(parisize, 4000000000);
print("=== (A) Atkin-Lehner data ===");
M4 = mfinit([7,4],0);
print("  dim S_4(Gamma_0(7))    = ", mfdim(mfinit([7,4],1)));
Mnew = mfinit([7,4],4);
print("  dim S_4^new(Gamma_0(7))= ", mfdim(Mnew));
F = mfeigenbasis(Mnew)[1];
print("  newform coefs 0..12    = ", mfcoefs(F,12));
{ iferr(print("  W_7 eigenvalue (mfatkineigenvalues) = ", mfatkineigenvalues(Mnew,7)), e, print("  atkin error: ", e)); }
print("");
print("  dim S_{5/2}(Gamma_0(28))            = ", mfdim(mfinit([28,5/2],1)));
K = mfkohnenbasis(mfinit([28,5/2],1));
print("  dim Kohnen plus subspace            = ", #K);
print("");
print("=== (B) T_{p^2} on the twisted-trace coefficients c(d), row s7 ===");
CDL = read("70_cd_s7.txt");
NCD = 400;
CD = vector(NCD);
KNOWN = vector(NCD);
{ for(i=1,#CDL, my(e=CDL[i]);
    if(e[1]<=NCD && type(e[2])!="t_STR", CD[e[1]]=e[2]; KNOWN[e[1]]=1)); }
{ adm(d) = if(d<=0, 0, if(d%4!=0 && d%4!=1, 0, my(s=((-3*d)%7+7)%7); s==0||s==1||s==2||s==4)); }
{ cc(d) = if(d<1 || d>NCD, 0, CD[d]); }
{ kn(d) = if(d<1 || d>NCD, 0, KNOWN[d]); }
\\ Weight 5/2 (k=2) Hecke on the index Delta = -3d :
\\   (f|T_{p^2})(Delta) = C(p^2 Delta) + (Delta/p) p C(Delta) + p^3 C(Delta/p^2)
\\ i.e. in terms of d:  c(p^2 d) + (-3d/p) p c(d) + p^3 c(d/p^2).
\\ TEST : is  c(p^2 d) + (-3d/p) p c(d) + p^3 c(d/p^2) - lam*c(d)  divisible by p^3,
\\ for lam = p (the value the theorem of 2.7 needs) and for lam = kronecker(-3,p)*p?
print("  d ranges over the admissible d with p^2 d <= ", NCD, "; c(49) is degenerate and skipped");
{ forprime(p=3,19,
   for(il=1,2,
     my(lam = if(il==1, p, kronecker(-3,p)*p), bad=0, cnt=0, first=0, worst=99);
     for(d=1, NCD\p^2,
       if(!adm(d), next);
       if(!kn(d) || !kn(p^2*d), next);
       if(d==49 || p^2*d==49, next);
       my(e = cc(p^2*d) + kronecker(-3*d,p)*p*cc(d) + if(d%(p^2)==0, p^3*cc(d/p^2), 0) - lam*cc(d));
       cnt++;
       my(v = if(e==0, 99, valuation(e,p)));
       if(v < 3, bad++; if(first==0, first=d));
       if(v < worst, worst = v));
     if(cnt>0,
       print("   p=",p,"  lam=", if(il==1,"p",Str(kronecker(-3,p),"*p")), "  tested ", cnt,
             " : ", if(bad==0, Str("p^3 | residual  PASS (min v_p = ", worst, ")"),
                                Str("FAIL ", bad, "/", cnt, " first d=", first, ", min v_p = ", worst))))));
}
print("");
print("=== control: same test on the level-one f4a coefficients (scalar, known-good) ===");
read("33_coeffs.gp");
{ ac(n) = if(n < -4 || n > NMAXC, 0, A4A[n+5]); }
{ forprime(p=3,13,
   my(bad=0,cnt=0,worst=99, lam=kronecker(-3,p)*p);
   for(n=1, NMAXC\p^2,
     my(e = ac(p^2*n) + kronecker(n,p)*p*ac(n) + if(n%(p^2)==0, p^3*ac(n/p^2), 0) - lam*ac(n));
     cnt++;
     my(v=if(e==0,99,valuation(e,p)));
     if(v<3, bad++); if(v<worst, worst=v));
   print("   p=",p," lam=chi_-3(p)*p : ", if(bad==0, Str("PASS (min v_p = ",worst,")"), Str("FAIL ",bad,"/",cnt))));
}
quit;
