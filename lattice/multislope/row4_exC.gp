default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
NMAX = #An - 1;
RR = 47; DD = 8; NCOL = (RR+1)*(DD+1);
kervec(pp) = {
  my(nlo = RR);
  my(nhi = min(NMAX, nlo+(NCOL\2)+40));
  my(nr = nhi-nlo+1);
  my(m1 = matrix(nr, NCOL, a, b, Mod(An[nlo+a-1-((b-1)\(DD+1))+1],pp) * Mod(nlo+a-1,pp)^((b-1)%(DD+1))));
  my(m2 = matrix(nr, NCOL, a, b, (Mod(numerator(Bn[nlo+a-1-((b-1)\(DD+1))+1]),pp)/Mod(denominator(Bn[nlo+a-1-((b-1)\(DD+1))+1]),pp)) * Mod(nlo+a-1,pp)^((b-1)%(DD+1))));
  my(kk = matker(matconcat([m1;m2])));
  if(matsize(kk)[2] != 1, return(0));
  my(v = kk[,1]);
  my(f = 0);
  for(t=1,NCOL, if(v[t]!=0, f=t; break));
  [f, v/v[f]]
}
{
  my(pl=List(), vl=List(), f0=0, pr=2^61-1, NP=8);
  for(t=1,NP,
    pr = precprime(pr-1);
    my(res = kervec(pr));
    if(res==0, print("kernel dim != 1"); quit);
    if(f0==0, f0=res[1], if(res[1]!=f0, print("pivot mismatch"); quit));
    listput(pl,pr); listput(vl,res[2])
  );
  print("pivot index = ", f0, "  (i=", (f0-1)\(DD+1), ", j=", (f0-1)%(DD+1), ")");
  cr = vector(NCOL);
  for(b=1,NCOL,
    my(cm = chinese(vector(NP, t, Mod(lift(vl[t][b]), pl[t]))));
    cr[b] = bestappr(cm)
  );
}
dd = 1; for(b=1,NCOL, dd = lcm(dd, denominator(cr[b])));
ciC = vector(NCOL, b, cr[b]*dd);
gg = 0; for(b=1,NCOL, gg = gcd(gg, ciC[b])); ciC = vector(NCOL, b, ciC[b]/gg);
QC = vector(RR+1, k, my(i=k-1); sum(j=0,DD, ciC[i*(DD+1)+j+1]*'n^j));
resid(x,n) = { my(sm=0); for(i=0,RR, if(n-i>=0, sm += subst(QC[i+1],'n,n)*x[n-i+1])); sm }
{my(bA=0,bB=0); for(n=0,NMAX, if(resid(An,n)!=0,bA++); if(resid(Bn,n)!=0,bB++));
 print("EXACT VERIFY n=0..",NMAX,": A failures=",bA,"  B failures=",bB);}
print("Q_0(n) = ", QC[1]);
print("factor Q_0 = ", factor(QC[1]));
print("Q_47(n) = ", QC[48]);
print("max coeff digits = ", vecmax(vector(NCOL,b,#Str(abs(ciC[b])))));
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recC.txt", ciC);
chC = sum(i=0,RR, ciC[i*(DD+1)+DD+1]*x^(RR-i));
chC = chC/content(chC);
print("char poly degree ", poldegree(chC));
print("chC factored: "); print(factor(chC));
print("product of char roots = ", (-1)^poldegree(chC)*polcoeff(chC,0)/pollead(chC));
print("factored product = ", factor((-1)^poldegree(chC)*polcoeff(chC,0)/pollead(chC)));
{for(pi=1,4, my(p=[2,3,5,7][pi]); print("newtonpoly(chC, p=",p,") = ", Vec(Set(newtonpoly(chC,p))), "  multiplicities: ", vector(#Set(newtonpoly(chC,p)), t, my(s=Vec(Set(newtonpoly(chC,p)))[t]); sum(u=1,poldegree(chC), if(newtonpoly(chC,p)[u]==s,1,0)))));}
quit
