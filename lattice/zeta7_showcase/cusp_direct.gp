default(parisize,4000000000);
default(realprecision,150);
N=5200;
ev(S,q0)=subst(truncate(S),q,q0);
d12=[1,2,3,4,6,12]; c12=[1,-572,11583,-36608,46332,-20736];
g=vector(N,m,sum(i=1,6,if(m%d12[i]==0,c12[i]*sigma(m\d12[i],7),0)));
PhiS=sum(m=1,N,g[m]*q^m)+O(q^(N+1));
PsiS=sum(m=1,N,(g[m]/m^7)*q^m)+O(q^(N+1));
xi=(209/1728)*zeta(7);
et(k)=eta(q^k+O(q^(N+2)));
hS=et(1)^3*et(4)*et(6)^2/(et(2)^2*et(3)*et(12)^3)/q;
xS=hS/((hS+3)*(hS+4));
DxS=q*deriv(xS,q);
pts=[0.02,0.017,0.015,0.013,0.011,0.010,0.009,0.008,0.007,0.0065,0.006,0.0055,0.005];
Ls=vector(#pts); Ws=vector(#pts);
{for(i=1,#pts,my(y=pts[i],tau=1/2+I*y,q0=exp(2*Pi*I*tau));
  my(x0=ev(xS,q0),Dx0=ev(DxS,q0),Phi0=ev(PhiS,q0),Psi0=ev(PsiS,q0));
  my(U1=(1+x0)*Phi0/Dx0, W1=U1*(Psi0-xi), L=log(1+x0));
  Ls[i]=L; Ws[i]=W1;
  print("y=",y,"  |q|=",abs(q0)*1.0,"  log(1+x)=",L,"  W'=",W1));}
\\ least squares polynomial fit of degree 7 in L
{my(m=#pts,Mt=matrix(m,8,i,j,Ls[i]^(j-1)),v=Ws~); my(sol=matsolve(Mt~*Mt,Mt~*v));
 print("alpha_7 = ",sol[8]);
 print("14161/5040 = ",14161/5040*1.0);
 print("alpha_6 = ",sol[7]); print("alpha_5 = ",sol[6]); print("alpha_0 = ",sol[1]);
 print("fit residual = ",sqrt(norml2(Mt*sol-v)/m));}
quit;
