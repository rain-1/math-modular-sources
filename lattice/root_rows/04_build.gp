\\ Build the A-rows of every Sym^w system with w>=2 in the project.
\\ Writes lattice/root_rows/rows_<name>.txt : one integer A_n per line.
default(seriesprecision, 400);
M = 520;
NN = 400;

ee(d)={my(r=1+O(q^(M+1))); for(n=1,(M\d)+1, r*=(1-q^(d*n))); r};

dump(name, A, n)={system(Str("rm -f lattice/root_rows/rows_",name,".txt"));
  for(i=1,n, write(Str("lattice/root_rows/rows_",name,".txt"), A[i]));
  print(name,": A_0..A_8 = ",vector(9,j,A[j]))};

q='q;
E1=ee(1);E2=ee(2);E3=ee(3);E4=ee(4);E6=ee(6);E8=ee(8);E12=ee(12);E16=ee(16);E24=ee(24);

\\ ---------- 1. zeta(7) level 24, w=6 :  A(z)=sqrt(1-34s+s^2)*Apery(s)^3, s=(z/(1-z))^2
{
 my(Ms=NN\2+4);
 my(ap=sum(m=0,Ms, sum(k=0,m,binomial(m,k)^2*binomial(m+k,k)^2)*'s^m)+O('s^(Ms+1)));
 my(Ac= ap^3*sqrt(1-34*'s+'s^2+O('s^(Ms+1))));
 my(zs='z+O('z^(NN+2)), ss=(zs/(1-zs))^2);
 my(Az=subst(Ac,'s,ss));
 my(A=vector(NN+1,i,polcoeff(Az,i-1)));
 dump("zeta7_L24", A, NN+1);
}

\\ ---------- 2. beta(4) level 24, w=3 :  A = F^3, F=eta2^10/(eta1^4 eta4^4)=theta3(q)^2,
\\              coordinate z = r/(1+r), r=(eta2 eta12/(eta4 eta6))^6
{
 my(F = E2^10/(E1^4*E4^4));
 my(th = sum(n=-20,20,q^(n^2))+O(q^(M+1)));
 print("beta4: F == theta3^2 ? ", F-th^2 == O(q^(M-2)));
 my(r = q*(E2*E12/(E4*E6))^6, z = r/(1+r));
 print("beta4: z = ",z+O(q^8));
 my(Qz = serreverse(z));
 my(Az = subst(F^3, q, Qz));
 my(A=vector(NN+1,i,polcoeff(Az,i-1)));
 dump("beta4_L24", A, NN+1);
 \\ also the weight-1 root directly in z, for cross-check
 my(gz = subst(F,q,Qz));
 system("rm -f lattice/root_rows/rows_beta4_L24_g.txt");
 for(i=1,NN+1, write("lattice/root_rows/rows_beta4_L24_g.txt", polcoeff(gz,i-1)));
 print("beta4: g_0..g_8 (direct) = ",vector(9,j,polcoeff(gz,j-1)));
}

\\ ---------- 3. L(4,chi_-3) level 24, w=3 : A = U(2tau), U=eta1^9/eta3^3=b(q)^3,
\\              coordinate x = eta1^2 eta24^2 eta4 eta6 /(eta2 eta12 eta3^2 eta8^2)
{
 my(U2 = E2^9/E6^3);              \\ = U(2tau) = b(q^2)^3
 my(b2 = E2^3/E6);                \\ = b(q^2)
 print("Lchi3: U(2t) == b(q^2)^3 ? ", U2-b2^3 == O(q^(M-2)));
 my(x = q*E1^2*E24^2*E4*E6/(E2*E12*E3^2*E8^2));
 print("Lchi3: x = ",x+O(q^8));
 my(Qx = serreverse(x));
 my(Ax = subst(U2,q,Qx));
 my(A=vector(NN+1,i,polcoeff(Ax,i-1)));
 dump("Lchi3_L24", A, NN+1);
 my(gx = subst(b2,q,Qx));
 system("rm -f lattice/root_rows/rows_Lchi3_L24_g.txt");
 for(i=1,NN+1, write("lattice/root_rows/rows_Lchi3_L24_g.txt", polcoeff(gx,i-1)));
 print("Lchi3: g_0..g_8 (direct) = ",vector(9,j,polcoeff(gx,j-1)));
}

\\ ---------- 4. zeta(5) level 16, w=4 : F = Phi16/Dt,  t = x/(8x^2+2x+1)
{
 my(X = q*E2*E16^2/(E1^2*E8), T = X/(8*X^2+2*X+1));
 my(cv=[1,-85,1428,-5440,4096], dv=[1,2,4,8,16]);
 my(g=vector(M)); for(n=1,M, my(s=0); for(i=1,5, if(n%dv[i]==0, s+=cv[i]*sigma(n\dv[i],5))); g[n]=s);
 my(PHI=sum(n=1,M,g[n]*q^n)+O(q^(M+1)));
 my(DT=q*deriv(T,q), F=PHI/DT);
 my(QT=serreverse(T), Az=subst(F,q,QT));
 my(A=vector(NN+1,i,polcoeff(Az,i-1)));
 dump("zeta5_L16", A, NN+1);
}

\\ ---------- 5. zeta(5) level 12, w=4 : two sources, Domb coordinate t=w/(1+w)^2
{
 my(W12=q*(E1*E12/(E3*E4))^4, T12=W12/(1+W12)^2);
 my(DT=q*deriv(T12,q), QT=serreverse(T12));
 my(dv=[1,2,3,4,6,12]);
 my(src=[[1,-104,351,832,-2808,1728],[1,-176,2079,-4928,4752,-1728]]);
 my(nms=["zeta5_L12p","zeta5_L12m"]);
 for(k=1,2,
   my(cv=src[k], g=vector(M));
   for(n=1,M, my(s=0); for(i=1,6, if(n%dv[i]==0, s+=cv[i]*sigma(n\dv[i],5))); g[n]=s);
   my(PHI=sum(n=1,M,g[n]*q^n)+O(q^(M+1)), F=PHI/DT);
   my(Az=subst(F,q,QT));
   my(A=vector(NN+1,i,polcoeff(Az,i-1)));
   dump(nms[k], A, NN+1));
}
\q
