\\ 63_table1.gp -- the master congruence n^2 | beta(n) across Pasol-Zudilin's Table 1 of
\\ strong magnetic weight-4 forms Psi(f_m), with psi = chi_{-m} (m = the pole discriminant).
M = 150;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
js  = E4s^3/Dl;
{ run(nm, D, sc, num, den) = my(Ph, c, cp, b, bad, s);
  Ph = sc*E4s*num/den^2;
  c = vector(M,m,polcoeff(Ph,m));
  bad = 0; for(m=1,M, if(denominator(c[m])!=1 || c[m]%m!=0, bad=m; break));
  print(nm,"  psi=chi_",D,"  c(1..3)=",vector(3,i,c[i]));
  print("    m | c(m), m<=",M,"? ", if(bad==0,"YES",concat("NO at ",bad)));
  if(bad!=0, return(0));
  cp = vector(M,m,c[m]/m);
  b = vector(M);
  for(n=1,M, my(t=0); fordiv(n,d, t += moebius(d)*kronecker(D,d)*cp[n/d]); b[n]=t);
  bad = List(); for(n=1,M, if(b[n]%(n^2)!=0, listput(bad,n)));
  print("    n^2 | beta(n)? ", if(#bad==0,"YES", concat(concat("NO, ",#bad)," failures: ")));
  if(#bad>0, print("       ",vector(min(12,#bad),i,Vec(bad)[i])));
  s=""; forprime(p=2,19, my(mn=99); for(n=1,M, if(n%p==0, my(v=valuation(b[n],p)-2*valuation(n,p)); if(v<mn,mn=v))); s=concat(s,concat(concat(" ",p),concat(":",mn))));
  print("    min v_p(beta)-2v_p(n):",s,"    beta(1..5)=",vector(5,i,b[i]));
}
run("Psi(f_7)  ",  -7,   27, 19*js-8*3375,          js+3375);
run("Psi(f_8)  ",  -8,   -8, 101*js-3*8000,         js-8000);
run("Psi(f_11) ", -11,   64, 43*js-6*32768,         js+32768);
run("Psi(f_19) ", -19, 1728, 25*js-2*884736,        js+884736);
run("Psi(f_43) ", -43, 1728, 11329*js-578*884736000, js+884736000);
run("Psi(f_15) ", -15,   15, 785*js^3-15219684*js^2+28709816985*js+837864*121287375, js^2+191025*js-121287375);
run("Psi(f_3|T4)",  -3,  48, 14*js+18*3375,         js-54000);
run("Psi(f_4|..)",  -4, 108, 611*js+404*35937,      js-287496);
quit;
