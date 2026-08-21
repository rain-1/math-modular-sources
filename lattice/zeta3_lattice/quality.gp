default(realprecision,40);
L2=log(2); LL=log(12+8*sqrt(2));   \\ log lambda_1^T
\\ Domb sampled at a*n, T sampled at c*n. k=3 denominator exponent.
\\ gainflag=1 : use 2-adic determinant divisor; 0 : archimedean only
Q(a,c,gainflag)={
 my(mx=max(a,c), S=3*mx, G=if(gainflag,min(6*a,4*c)*L2,0));
 my(sig=S+G);
 my(A1=S+4*a*L2, E1=S+2*a*L2, A2=S+c*LL, E2=S+4*c*L2-c*LL);
 my(gap=(A2+E1)-(A1+E2));
 my(x=(sig+E2-E1)/2);
 my(H=A2-x, F=x+E1-sig);
 [H, F, (H-F)/H, gap, sig];
};
print("== quality delta(a,c) with 2-adic gain ==");
{my(best=0,ba,bc);
 for(a=1,12,for(c=1,18, my(r=Q(a,c,1)); if(r[4]>0 && r[1]>0 && r[3]>best, best=r[3];ba=a;bc=c)));
 print("best over integer a<=12,c<=18: a:c=",ba,":",bc,"  delta=",best);}
{my(best=0,br); forstep(r=0.4,4,0.001, my(q=Q(1,r,1)); if(q[4]>0&&q[1]>0&&q[3]>best,best=q[3];br=r));
 print("continuous optimum r=c/a=",br,"  delta=",best);}
print();
{my(r=Q(2,3,1)); print("a:c = 2:3 (the 2-adic balance point):");
 print("  sigma=",r[5],"  H=",r[1],"  F=",r[2],"  gap=",r[4]);
 print("  delta=",r[3]);
 print("  irrationality measure mu <= H/(H-F) = ",r[1]/(r[1]-r[2]));}
print();
print("== archimedean only (no hidden 2-adic determinant) ==");
{my(best=0,ba,bc);for(a=1,12,for(c=1,18,my(r=Q(a,c,0));if(r[4]>0&&r[1]>0&&r[3]>best,best=r[3];ba=a;bc=c)));
 print("best a:c=",ba,":",bc,"  delta=",best);}
{my(r=Q(2,3,0)); print("a:c=2:3 arch-only: H=",r[1]," F=",r[2]," delta=",r[3]);}
{my(best=0,br); forstep(r=0.4,4,0.001, my(q=Q(1,r,0)); if(q[4]>0&&q[1]>0&&q[3]>best,best=q[3];br=r));
 print("continuous arch-only optimum r=",br," delta=",best);}
print();
print("== symbolic check at 2:3 ==");
print("  H = (9 + 9*logLambda - 20*log2)/2 = ",(9+9*LL-20*L2)/2);
print("  F = (9 + 20log2 - 3logLambda)/2 - 8log2 = ",(9+20*L2-3*LL)/2-8*L2);
print("  H-F = 6logLambda-12log2 = ",6*LL-12*L2);
print("  delta = (12logLambda-24log2)/(9+9logLambda-20log2) = ",(12*LL-24*L2)/(9+9*LL-20*L2));
print();
print("== sensitivity: what denominator exponent k would be needed for delta>=1 at 2:3? ==");
{Qk(a,c,k)=my(mx=max(a,c),S=k*mx,G=min(6*a,4*c)*L2,sig=S+G,A1=S+4*a*L2,E1=S+2*a*L2,A2=S+c*LL,E2=S+4*c*L2-c*LL,x=(sig+E2-E1)/2,H=A2-x,F=x+E1-sig);(H-F)/H;
 forstep(k=0,3,0.05, print1("k=",k," d=",Qk(2,3,k),"  ")); print();}
\q
