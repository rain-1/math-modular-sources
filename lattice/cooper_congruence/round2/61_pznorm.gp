\\ 61_pznorm.gp -- with Pasol-Zudilin's own normalisations (64 F_4a, 108 F_4b: the constants
\\ that make the weight-5/2 input integral) AND the genus character of the pole discriminant,
\\ does the master congruence n^2 | beta(n) hold at EVERY prime?
M = 500;
E4s = 1 + 240*sum(n=1,M,sigma(n,3)*'q^n) + O('q^(M+1));
E6s = 1 - 504*sum(n=1,M,sigma(n,5)*'q^n) + O('q^(M+1));
Dl  = 'q*prod(n=1,M,(1-'q^n+O('q^(M+1)))^24);
{ betav(S,D,sc) = my(c=vector(M,m,sc*polcoeff(S,m)), cp, b=vector(M));
  cp = vector(M,m,c[m]/m);
  for(n=1,M, my(s=0); fordiv(n,d, s += moebius(d)*if(D==0,1,kronecker(D,d))*cp[n/d]); b[n]=s);
  b; }
{ rep(nm,S,D,sc) = my(b=betav(S,D,sc), bad=List());
  for(n=1,M, if(b[n]%(n^2)!=0, listput(bad,n)));
  print(nm,"  psi=chi_",D,"  scale=",sc,":   n^2 | beta(n) for all n<=",M,"?  ",
        if(#bad==0,"YES", concat(concat("NO, ",#bad)," failures, first: ")));
  if(#bad>0, print("     ", vector(min(15,#bad),i,Vec(bad)[i])));
  print("     beta(1..6)=",vector(6,i,b[i]));
  print("     min over n<=",M," of v_p(beta(n))-2v_p(n), p<=29:");
  my(s="");
  forprime(p=2,29, my(mn=999); for(n=1,M, if(n%p==0, my(v=valuation(b[n],p)-2*valuation(n,p)); if(v<mn,mn=v))); s=concat(s,concat(concat(" ",p),concat(":",mn))));
  print("    ",s);
}
rep("F_4a", Dl/E4s^2, -3, 64);
rep("F_4b", E4s*Dl/E6s^2, -4, 108);
rep("F_4a", Dl/E4s^2, -3, 1);
rep("F_4b", E4s*Dl/E6s^2, -4, 1);
quit;
