default(parisize, 8000000000);
default(realprecision, 200);
OUT="/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/";
{P3 = [nn^2+2*nn+1, -5*nn^2-5*nn-2, -4*nn^2+16*nn-4, 40*nn^2-120*nn+76, -16*nn^2+32*nn+32, -80*nn^2+560*nn-992, 64*nn^2-512*nn+1024];}
{P6 = [nn^2+2*nn+1, 9*nn^2+nn, 20*nn^2-40*nn+20, -20*nn^2+20*nn+4, -96*nn^2+352*nn-336, -64*nn^2+320*nn-400];}
{ext(QQ, init, NMX) =   /* init = vector of u_0..u_{r} (length r+1) ; extend */
  my(r=#QQ-1, U=vector(NMX+1));
  for(i=1,#init, U[i]=init[i]);
  for(n=#init-1, NMX-1, my(s=0);
    for(j=1,r, my(idx=n+1-j); if(idx>=0, s += subst(QQ[j+1],nn,n)*U[idx+1]));
    U[n+2] = -s/subst(QQ[1],nn,n));
  U;}
{canon(QQ, NMX) = my(r=#QQ-1, U=vector(NMX+1)); U[2]=1;
  for(n=1,NMX-1, my(s=0);
    for(j=1,r, my(idx=n+1-j); if(idx>=0, s += subst(QQ[j+1],nn,n)*U[idx+1]));
    U[n+2] = -s/subst(QQ[1],nn,n));
  U;}
{doit(QQ, afile, tag, NMX) =
  my(r=#QQ-1, Af=readvec(afile), NA=#Af-1);
  my(A = ext(QQ, vector(r+1,i,Af[i]), NMX));
  print("\n##### ", tag);
  print("  a_n check vs file to n=",min(NA,NMX),": ", vector(min(NA,NMX)+1,i,A[i])==vector(min(NA,NMX)+1,i,Af[i]));
  my(W = canon(QQ, NMX));
  print("  w_n first 8: ", vector(8,i,W[i]));
  foreach([50,100,200,300,400,600,800,NMX], m, if(m<=NMX,
     print("   xi(",m,") = ", W[m+1]/A[m+1]*1.0)));
  print("   |xi(NMX)-xi(NMX-1)| = ", abs(W[NMX+1]/A[NMX+1]*1.0 - W[NMX]/A[NMX]*1.0));
  W[NMX+1]/A[NMX+1]*1.0;}
x3 = doit(P3, Str(OUT,"lvl12_p3_A.txt"), "placement 3 (c=3)", 800);
x6 = doit(P6, Str(OUT,"lvl12_p6_A.txt"), "placement 6 (c=6)", 800);
print("\n-G/27 = ", -Catalan/27, "   -(4/9)z2 = ", -4*zeta(2)/9);
print("G/9 = ", Catalan/9, "   (2/9)z2 = ", 2*zeta(2)/9);
FN=Str(OUT,"lvl12_limits.txt"); write1(FN,"");
write(FN, Str("xi_canon_lvl12_p3 ", x3));
quit;
