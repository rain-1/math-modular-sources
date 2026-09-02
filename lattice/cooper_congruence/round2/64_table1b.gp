\\ 64_table1b.gp -- the rest of Pasol-Zudilin's Table 1.
M = 120;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
js  = E4s^3/Dl;
{ run(nm, D, sc, num, den) = my(Ph, c, cp, b, bad, s);
  Ph = sc*E4s*num/den^2;
  c = vector(M,m,polcoeff(Ph,m));
  bad = 0; for(m=1,M, if(denominator(c[m])!=1 || c[m]%m!=0, bad=m; break));
  print(nm,"  psi=chi_",D,"   m|c(m)? ", if(bad==0,"YES",concat("NO at ",bad)));
  if(bad!=0, return(0));
  cp = vector(M,m,c[m]/m);
  b = vector(M);
  for(n=1,M, my(t=0); fordiv(n,d, t += moebius(d)*kronecker(D,d)*cp[n/d]); b[n]=t);
  bad = List(); for(n=1,M, if(b[n]%(n^2)!=0, listput(bad,n)));
  print("    n^2 | beta(n), n<=",M,"? ", if(#bad==0,"YES", concat(concat("NO, ",#bad)," failures: ")));
  if(#bad>0, print("       ",vector(min(12,#bad),i,Vec(bad)[i])));
  s=""; forprime(p=2,23, my(mn=99); for(n=1,M, if(n%p==0, my(v=valuation(b[n],p)-2*valuation(n,p)); if(v<mn,mn=v))); s=concat(s,concat(concat(" ",p),concat(":",mn))));
  print("    min v_p(beta)-2v_p(n):",s);
}
run("Psi(f_67) ", -67, 1728, 1221961*js-49442*5280^3, js+5280^3);
run("Psi(f_163)",-163, 1728, 908855380249*js-23238932978*640320^3, js+640320^3);
run("Psi(f_20) ", -20,  -80, 733*js^3+72767680*js^2-984198615040*js+12^3*20^3*880^3, js^2-158*20^3*js-880^3);
run("Psi(f_7|.)",  -7,   27, 82451*js+5272*255^3, js-255^3);
{ my(P23);
  P23 = 141826*js^5 - 286458244*125*js^4 + 5214621227*15625*js^3 + 3414887843776*1953125*js^2 - 47816219216827*244140625*js + 4378632*1873*1953125*125;
  run("Psi(f_23) ", -23, -1, P23, js^3+27934*125*js^2-329683*15625*js+187^3*1953125);
}
quit;
