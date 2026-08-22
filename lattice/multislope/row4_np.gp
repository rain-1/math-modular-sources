default(parisizemax, 6000000000);
ciA = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recA.txt");
ciB = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt");
ciC = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recC.txt");
QC = vector(48, k, my(i=k-1); sum(j=0,8, ciC[i*9+j+1]*'n^j));
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_QC.txt", QC);
mkch(ci,r,D) = sum(i=0,r, ci[i*(D+1)+D+1]*x^(r-i));
chA = mkch(ciA,22,23); chA=chA/content(chA);
chB = mkch(ciB,35,7);  chB=chB/content(chB);
chC = mkch(ciC,47,8);  chC=chC/content(chC);
q5  = (x-1)*(x^2-8*x+8)*(x^2+4*x-4);
npsum(f,p) = { my(v=newtonpoly(f,p)); my(s=Set(v)); vector(#s, t, [s[t], sum(u=1,#v, if(v[u]==s[t],1,0))]) }
{for(nm=1,4,
  my(f=[q5,chA,chB,chC][nm]);
  print("--- ", ["base quintic (x-1)(x^2-8x+8)(x^2+4x-4)","chA r=22 D=23","chB r=35 D=7","chC r=47 D=8"][nm],
        "  deg=",poldegree(f), "  prod roots=", (-1)^poldegree(f)*polcoeff(f,0)/pollead(f));
  for(pi=1,4, my(p=[2,3,5,7][pi]);
    print("   p=",p,"  [slope,mult] = ", npsum(f,p))));}
print("");
print("v_p of product of roots:");
{for(nm=1,4, my(f=[q5,chA,chB,chC][nm]); my(pr=(-1)^poldegree(f)*polcoeff(f,0)/pollead(f));
  print("  ", ["quintic","chA","chB","chC"][nm], ": ", pr, "  v_2=",valuation(pr,2)," v_3=",valuation(pr,3)," v_5=",valuation(pr,5)," v_7=",valuation(pr,7)));}
quit
