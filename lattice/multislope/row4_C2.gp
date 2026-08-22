default(parisizemax, 6000000000);
default(realprecision,60);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
ciC = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recC.txt");
RR=47; DD=8;
QC = vector(RR+1, k, my(i=k-1); sum(j=0,DD, ciC[i*(DD+1)+j+1]*'n^j));
resid(x,n) = { my(sm=0); for(i=0,RR, if(n-i>=0, sm += subst(QC[i+1],'n,n)*x[n-i+1])); sm }
print("n where L_C[B] != 0 : ", select(k->resid(Bn,k)!=0, vector(1598,k,k-1)));
print("extra char factor roots: ", polroots(952*x^7 - 9519*x^6 + 31599*x^5 - 41048*x^4 - 159840*x^3 + 877656*x^2 - 1071592*x + 322560));
print("is it irreducible? ", polisirreducible(952*x^7 - 9519*x^6 + 31599*x^5 - 41048*x^4 - 159840*x^3 + 877656*x^2 - 1071592*x + 322560));
NC = 260;
companion(s) = { my(x=vector(NC+1)); x[s+1]=1; for(n=s+1,NC, my(sm=0); for(i=1,RR, if(n-i>=0, sm += subst(QC[i+1],'n,n)*x[n-i+1])); x[n+1] = -sm/(952*n^8)); x }
Y1=companion(1); Y2=companion(2);
print("Y1[0..6] = ", vector(7,i,Y1[i]));
print("Bn[0..6] = ", vector(7,i,Bn[i]));
{my(ok=1); for(n=0,NC, if(Y1[n+1]!=Bn[n+1], ok=0; print("  first differ at n=",n); break)); print("Y1==B_n for n<=",NC,"? ",ok);}
print("ratio Y1_n/B_n at n=10,20,30: ", vector(3,k,my(n=10*k); if(Bn[n+1]!=0,1.0*Y1[n+1]/Bn[n+1],"NA")));
sharpk(x,nm)={my(kb=0,dn=1); for(n=1,nm, dn=lcm(dn,n); my(k=0); while(denominator(dn^k*x[n+1])!=1, k++; if(k>40,break)); if(k>kb,kb=k)); kb}
print("Y^(1) sharp k = ", sharpk(Y1,150));
print("Y^(2) sharp k = ", sharpk(Y2,150));
{my(k=sharpk(Y1,150),dn=1,bad=0); for(n=1,150, dn=lcm(dn,n); if(denominator(dn^(k-1)*Y1[n+1])!=1,bad++)); print("  Y1: k-1=",k-1," fails at ",bad,"/150");}
{my(k=sharpk(Y2,150),dn=1,bad=0); for(n=1,150, dn=lcm(dn,n); if(denominator(dn^(k-1)*Y2[n+1])!=1,bad++)); print("  Y2: k-1=",k-1," fails at ",bad,"/150");}
{my(ns=[20,40,60,80,100,120,140,160,180,200,220,240,260]);
 for(w=1,2, my(x=if(w==1,Y1,Y2)); print("companion Y^(",w,")");
   for(pi=1,4, my(p=[2,3,5,7][pi]);
     print("  p=",p," v_p(x_n/A_n - x_{n-1}/A_{n-1}) at n=",ns," : ", vector(#ns,k,my(n=ns[k]); my(z=x[n+1]/An[n+1]-x[n]/An[n]); if(z==0,"inf",valuation(z,p))))));}
print("Y1_n/A_n at n=260 = ", 1.0*Y1[261]/An[261]);
print("Y2_n/A_n at n=260 = ", 1.0*Y2[261]/An[261]);
print("(1463/13824)zeta(7) = ", 1.0*(1463/13824)*zeta(7));
print("arch decay |Y1_n/A_n - Y1_{n-1}/A_{n-1}|^(1/n) at n=260: ", abs(1.0*(Y1[261]/An[261]-Y1[260]/An[260]))^(1.0/260));
print("arch decay for B: ", abs(1.0*(Bn[261]/An[261]-Bn[260]/An[260]))^(1.0/260));
quit
