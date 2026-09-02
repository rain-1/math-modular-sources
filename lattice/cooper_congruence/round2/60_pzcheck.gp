\\ 60_pzcheck.gp -- INDEPENDENT CHECK of the claim that Pasol-Zudilin's level-one magnetic
\\ forms satisfy the master congruence n^2 | beta(n) once the character is the genus
\\ character of the POLE discriminant: psi = chi_{-3} for F_4a (pole at rho), chi_{-4} for
\\ F_4b (pole at i).  Round 1 section 3 tested psi = 1 and concluded "refuted".
M = 400;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
E6s = 1 - 504*sum(n=1,M,sigma(n,5)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
F4a = Dl/E4s^2;
F4b = E4s*Dl/E6s^2;
F6  = E6s*Dl/E4s^3;
{ betav(S,D) = my(c=vector(M,m,polcoeff(S,m)), cp, b=vector(M));
  cp = vector(M,m,c[m]/m);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*if(D==0,1,kronecker(D,d))*cp[n/d]); b[n]=s);
  b; }
{ rep(nm, S, D) = my(b=betav(S,D), bad=List(), v);
  for(n=2,M, if(b[n]%(n^2)!=0, listput(bad,n)));
  print(nm, "  psi = ", if(D==0,"1",concat("chi_",D)), " :  #{n<=",M,": n^2 nmid beta(n)} = ", #bad);
  print("      first failures: ", vector(min(25,#bad),i,Vec(bad)[i]));
  print("      beta(1..8) = ", vector(8,i,b[i]));
}
rep("F_4a", F4a, 0);
rep("F_4a", F4a, -3);
rep("F_4b", F4b, 0);
rep("F_4b", F4b, -4);
\\ which primes are responsible?
{ pfail(nm,S,D) = my(b=betav(S,D));
  print(nm," psi=",D,"  v_p(beta(p^k)) - 2k  for p<=31, k=1,2:");
  forprime(p=2,31,
    my(s="");
    for(k=1,2, if(p^k<=M, s = concat(s, concat(concat(" ", valuation(b[p^k],p)-2*k),""))));
    print("   p=",p,"  ",s));
}
pfail("F_4a", F4a, -3);
pfail("F_4b", F4b, -4);
quit;
