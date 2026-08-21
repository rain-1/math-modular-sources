default(parisizemax,2000000000); default(realprecision,60);
M=160; default(seriesprecision,M+4);
et(k)=eta(q^k+O(q^(M+4)));
F=et(2)^7*et(3)^7/(et(1)^5*et(6)^5); t=q*(et(1)*et(6)/(et(2)*et(3)))^12;
f=sqrt(F); u=t/4; Dq(g)=q*deriv(g,q);
Phi=f*Dq(u);
print("Phi = ",vector(12,i,polcoeff(Phi,i)));
\\ is Phi = Phi_gamma/(4 f) ?  Phi_gamma = F*Dq(t)
print("check Phi*4*f == F*Dq(t): ", Phi*4*f-F*Dq(t)==0);
\\ denominators of Phi coefficients
print("denominators: ",vector(20,i,denominator(polcoeff(Phi,i))));
\\ try mftobasis in M_3(Gamma0(N),chi) for several N, chi
cands=[[24,-24],[24,-8],[24,-3],[24,-4],[48,-24],[48,-8],[48,-3],[48,-4],[12,-3],[12,-4],[72,-24],[72,-8],[72,-3]];
{for(i=1,#cands, my(Nl=cands[i][1],chi=cands[i][2]); my(mf=mfinit([Nl,3,chi],4)); my(d=mfdim(mf)); if(d==0,next);
 my(sb=mfsturm(mf)); my(vv=vector(sb+2,j,polcoeff(Phi,j-1))); my(ok=1);
 my(bb=mftobasis(mf,vv,1)); if(bb==0||type(bb)!="t_COL", print("N=",Nl," chi=",chi," dim ",d,": no"), print("N=",Nl," chi=",chi," dim ",d,": FITS (sturm ",sb,")")));}
\q
