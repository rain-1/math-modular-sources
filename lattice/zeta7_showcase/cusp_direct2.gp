default(realprecision,160);
d12=[1,2,3,4,6,12]; c12=[1,-572,11583,-36608,46332,-20736];
xi=(209/1728)*zeta(7);
eps=10.0^(-150);
etaq(Q)={my(s=0,m=0,t); while(1, t=Q^(m*(3*m-1)/2); if(abs(t)<eps,break); s+=(-1)^m*(t+if(m>0,Q^(m*(3*m+1)/2),0)); m++); s};
\\ eta(tau)=q^{1/24} prod(1-q^n) = q^{1/24} * etaq(q) ; we only need quotients with total q^{-1} power handled explicitly
E2(Q)={my(s=0,n=1,t); while(1, t=n*Q^n/(1-Q^n); if(abs(t)<eps,break); s+=t; n++); 1-24*s};
E8(Q)={my(s=0,n=1,t); while(1, t=n^7*Q^n/(1-Q^n); if(abs(t)<eps,break); s+=t; n++); 1+480*s};
Lam7(Q)={my(s=0,e=1,t); while(1, t=Q^e/(1-Q^e)/e^7; if(abs(t)<eps && e>50,break); s+=t; e++); s};
rd=[[1,3],[4,1],[6,2],[2,-2],[3,-1],[12,-3]];
{doit(y)=my(tau=1/2+I*y, q=exp(2*Pi*I*tau));
  my(h=q^(-1)*prod(i=1,#rd, etaq(q^rd[i][1])^rd[i][2]));   \\ h = eta_1^3 eta_4 eta_6^2/(eta_2^2 eta_3 eta_12^3) = q^{-1} prod(1-q^..)
  my(Dlogh=sum(i=1,#rd, rd[i][2]*rd[i][1]/24*E2(q^rd[i][1])));
  my(x=h/((h+3)*(h+4)));
  my(dxdh=((h+3)*(h+4)-h*(2*h+7))/((h+3)*(h+4))^2);
  my(Dx=dxdh*h*Dlogh);
  my(Phi=sum(i=1,6,c12[i]*E8(q^d12[i]))/480);
  my(Psi=sum(i=1,6,c12[i]/d12[i]^7*Lam7(q^d12[i])));
  my(U1=(1+x)*Phi/Dx, W1=U1*(Psi-xi), L=log(1+x));
  [L,W1,x]};
pts=[0.02,0.015,0.012,0.010,0.008,0.007,0.006,0.005,0.0045,0.004,0.0035,0.003,0.0025,0.002];
Ls=vector(#pts); Ws=vector(#pts);
{for(i=1,#pts, my(r=doit(pts[i])); Ls[i]=r[1]; Ws[i]=r[2]; print("y=",pts[i],"  L=",real(r[1]),"  imL=",imag(r[1]),"  W'/L^7=",real(r[2]/r[1]^7)));}
{my(m=#pts,Mt=matrix(m,8,i,j,Ls[i]^(j-1)),v=Ws~); my(sol=matsolve(Mt~*Mt,Mt~*v));
 print("alpha_7 = ",real(sol[8]),"   imag part ",imag(sol[8]));
 print("14161/5040 = ",14161/5040*1.0);
 print("alpha_6 = ",real(sol[7]),"  alpha_5 = ",real(sol[6]));
 print("fit residual/|W| = ",sqrt(norml2(Mt*sol-v)/m)/abs(Ws[#pts]));}
quit;
