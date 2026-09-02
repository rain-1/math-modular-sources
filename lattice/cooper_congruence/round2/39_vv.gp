\\ 39_vv.gp -- the vector-valued reading of Cooper's row s7 at level 28.
\\ Uses the twisted CM trace coefficients c(d) (companion strand, 70_cd_s7.txt),
\\ supported on the admissible set  -3d = square mod 28
\\ (equivalently d = 0,1 mod 4 AND d = 0,1,2,4 mod 7).
\\
\\ (A) the Fricke/Atkin-Lehner obstruction at level 7;
\\ (B) the divisibility  p^{2a} | c(p^{2a} d)  for EVERY admissible d, not only squares;
\\ (C) the Hecke fit  c(p^2 d) + tau*p*c(d) + p^3 c(d/p^2) = lambda c(d)  mod p^e, tau=+1.
default(parisize, 4000000000);
print("================ (A) Fricke data at level 7 ================");
print("  dim S_4(Gamma_0(7)) = ", mfdim(mfinit([7,4],1)));
Mnew = mfinit([7,4],4);
F = mfeigenbasis(Mnew)[1];
CF = mfcoefs(F,12);
print("  weight-4 newform on Gamma_0(7), coefs 0..12 = ", CF);
print("  a_7 = ", CF[8], " => classical Fricke eigenvalue eps = -a_7/7^(k/2-1) = ", -CF[8]/7);
print("  (PARI mfatkineigenvalues = ", mfatkineigenvalues(Mnew,7), ")");
print("  Cooper's source has Phi|_4 W_7 = -Phi, so it lives in the eps = -1 eigenspace,");
print("  and dim S_4(Gamma_0(7))^{W_7 = -1} = 0.");
print("  ==> the uniqueness input of PZ's Lemma 2 (principal part determines the form)");
print("      IS available at level 28 on the W_7-antiinvariant part.");
CDL = read("70_cd_s7.txt");
NCD = 400;
CD = vector(NCD);
KN = vector(NCD);
{ for(i=1,#CDL, my(e=CDL[i]); if(e[1]<=NCD && type(e[2])!="t_STR", CD[e[1]]=e[2]; KN[e[1]]=1)); }
{ adm(d) = if(d<=0,0,if(d%4!=0 && d%4!=1,0, my(s=((-3*d)%7+7)%7); s==0||s==1||s==2||s==4)); }
{ cc(d) = if(d<1||d>NCD,0,CD[d]); }
print("");
print("================ (B)  p^{2a} | c(p^{2a} d)  for ALL admissible d ================");
print("  (extends  n^2 | beta_s7(n) = c(n^2)  from square indices to every index)");
{ forprime(p=2,19,
   for(a=1,3,
     my(bad=0, cnt=0, mn=99, lst=[]);
     if(p^(2*a) > NCD, break);
     for(d=1, NCD\p^(2*a),
       if(!adm(d) || d%p==0, next);
       if(!KN[d] || !KN[p^(2*a)*d], next);
       if(d==49 || p^(2*a)*d==49, next);
       cnt++;
       my(x = cc(p^(2*a)*d), v = if(x==0, 999, valuation(x,p)));
       if(v < 2*a, bad++; lst=concat(lst,[d]));
       if(v-2*a < mn, mn = v-2*a));
     if(cnt>0, print("   p=",p," a=",a," tested ",cnt," admissible d : ",
        if(bad==0, Str("PASS  (min excess ", mn, ")"), Str("FAIL ", bad, "/", cnt, " at d=", lst))))));
}
print("");
print("================ (C) Hecke fit with middle coefficient tau = +1 ================");
{ forprime(p=2,19,
   for(e=2,3,
     my(vals=[]);
     for(d=1, NCD\p^2,
       if(!adm(d) || !KN[d] || !KN[p^2*d] || d==49 || p^2*d==49, next);
       if(CD[d]==0 || CD[d]%p==0, next);
       my(lhs = cc(p^2*d) + p*cc(d) + if(d%(p^2)==0, p^3*cc(d/p^2), 0));
       vals = concat(vals, [lift(Mod(lhs,p^e)/Mod(CD[d],p^e))]));
     if(#vals>0,
       print("   p=",p," mod p^",e," (",#vals," values) : ",
             if(#Set(vals)==1, Str("CONSTANT lambda = ", vals[1], if(vals[1]==p%p^e, "  =  p   ***", "")), Str("not constant: ", vals))))));
}
print("");
print("  tau = +1 means the middle Hecke term carries the TRIVIAL character: the genus");
print("  character chi_{-3} of the twist has squared away, chi_{-3}(p)^2 = 1.  That is");
print("  exactly the mechanism predicted in FINDINGS_PZ.md section 3.");
print("  CAVEAT: every tested d is coprime to p, so the vanishing of the middle term when");
print("  p | d -- which a genuine T_{p^2} requires -- is NOT tested here.");
quit;
