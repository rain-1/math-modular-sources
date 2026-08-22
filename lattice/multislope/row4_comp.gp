default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
ciB = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt");
rB = 35; DB = 7;
QB = vector(rB+1, k, my(i=k-1); sum(j=0,DB, ciB[i*(DB+1)+j+1]*'n^j));
print("Q_0(n) = ", QB[1]);
print("Q_1(n) = ", QB[2]);

NC = 260;
companion(s) = {
  my(x = vector(NC+1));
  for(t=0, NC, x[t+1] = 0);
  x[s+1] = 1;
  for(n = s+1, NC,
    my(sm = 0);
    for(i = 1, rB, if(n-i >= 0, sm += subst(QB[i+1],'n,n)*x[n-i+1]));
    x[n+1] = -sm / n^7
  );
  x
}
X1 = companion(1);
X2 = companion(2);

print("\n--- check A_n satisfies the recurrence at every n>=0 (x_j=0 for j<0) ---");
{my(bad=0); for(n=0, 1000, my(sm=0); for(i=0,rB, if(n-i>=0, sm += subst(QB[i+1],'n,n)*An[n-i+1])); if(sm!=0, bad++));
 print("failures n=0..1000: ", bad);}

print("\n--- s=1 companion vs B_n ---");
print("X1[0..8] = ", vector(9,i,X1[i]));
print("B_n[0..8] = ", vector(9,i,Bn[i]));
{my(ok=1); for(n=0,NC, if(X1[n+1] != Bn[n+1], ok=0; if(n<40,print("  differ at n=",n,": X1=",X1[n+1]," B=",Bn[n+1])); break));
 print("X1 == B_n exactly for n<=",NC,"? ", ok);}
{my(rat = vector(6,i,my(n=i+9); if(Bn[n+1]!=0, X1[n+1]/Bn[n+1], "NA"))); print("ratios X1_n/B_n at n=10..15: ", rat);}

print("\n--- sharp lcm exponent k ---");
sharpk(x, nm) = {
  my(kb = 0);
  for(n=1, nm,
    my(dn = 1); for(t=2,n, dn=lcm(dn,t));
    my(d = denominator(x[n+1]));
    if(d > 1,
      my(k=0); while(denominator(dn^k * x[n+1]) != 1, k++; if(k>40, break));
      if(k > kb, kb = k)
    )
  );
  kb
}
print("X^(1): sharp k = ", sharpk(X1, 150));
print("X^(2): sharp k = ", sharpk(X2, 150));
{my(k=sharpk(X1,150)); my(dn=1,bad=0); for(n=1,150, dn=lcm(dn,n); if(denominator(dn^(k-1)*X1[n+1])!=1, bad++));
 print("  X^(1): k-1=",k-1," fails at ",bad," of 150 n"); }
{my(k=sharpk(X2,150)); my(dn=1,bad=0); for(n=1,150, dn=lcm(dn,n); if(denominator(dn^(k-1)*X2[n+1])!=1, bad++));
 print("  X^(2): k-1=",k-1," fails at ",bad," of 150 n"); }

print("\n--- p-adic valuations v_p(x_n/A_n - x_{n-1}/A_{n-1}) ---");
{my(ns=[20,40,60,80,100,120,140,160,180,200,220,240,260]);
 for(w=1,2,
   my(x = if(w==1, X1, X2));
   print("companion X^(",w,")");
   for(pi=1,4, my(p=[2,3,5,7][pi]);
     my(vv = vector(#ns, k, my(n=ns[k]); my(z = x[n+1]/An[n+1] - x[n]/An[n]); if(z==0, "inf", valuation(z,p))));
     print("  p=",p," at n=",ns," : ", vv)
   )
 );}
print("\n--- v_p(A_n) ---");
{my(ns=[20,60,100,140,180,220,260]);
 for(pi=1,4, my(p=[2,3,5,7][pi]); print("  v_",p,"(A_n) at n=",ns,": ", vector(#ns,k,valuation(An[ns[k]+1],p))));}

print("\n--- Newton polygons of char polys ---");
ch5 = (x-1)*(x^2-8*x+8)*(x^2+4*x-4);
print("base quintic = ", ch5, "  disc-free product of roots = ", -polcoeff(ch5,0)/pollead(ch5)*(-1)^0);
{for(pi=1,4, my(p=[2,3,5,7][pi]);
  print("  p=",p," newtonpoly(quintic) = ", newtonpoly(ch5,p)));}
chB = sum(i=0,rB, ciB[i*(DB+1)+DB+1]*x^(rB-i));
print("  chB == quintic^7 ? ", chB == ch5^7);
{for(pi=1,4, my(p=[2,3,5,7][pi]);
  print("  p=",p," newtonpoly(chB, deg 35) = ", newtonpoly(chB,p)));}
quit
