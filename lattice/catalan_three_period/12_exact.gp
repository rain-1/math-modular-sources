/* 12_exact.gp -- recognise the period polynomials (alp,bet) at the fold cusp
   x=-1/4 (the cusp 1/2 of Gamma_0(16), width 4, gamma_1=[-7,4;-16,9]) as exact
   Q-combinations of  G=Catalan (inner)  /  zeta(2) and i*pi^2 ;
   then solve the fold-regularity and target-zero conditions on both
   orientations of the 6-dimensional chi_{-4} weight-3 Eisenstein space.     */
default(parisizemax, 8*10^9);
default(realprecision, 130);
PRECDIG = 130;
NTOP = 4000;
chim(n) = if(n%2==0, 0, if(n%4==1, 1, -1));
ainn = vector(NTOP, n, sumdiv(n, d, chim(n/d)*d^2));
aout = vector(NTOP, n, sumdiv(n, d, chim(d)*d^2));
shft(v, dd) = vector(#v, n, if(n%dd==0, v[n/dd], 0));
thv(v, tau, NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n^2*p)); s;
thd(v, tau, NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n*p));  2*Pi*I*s;
nterms(tau) = my(y=imag(tau)); ceil((PRECDIG+25)*log(10)/(2*Pi*y));
gsh(gm,tau) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2],ch=chim(d)); 2*Pi^2*((a*tau+b)^2/(ch*(c*tau+d)) - tau^2);
gshd(gm,tau)= my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2],ch=chim(d)); 2*Pi^2*((2*a*(a*tau+b)*(c*tau+d)-c*(a*tau+b)^2)/(ch*(c*tau+d)^2)-2*tau);
{Rval(gm,v,tau) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2],ch=chim(d));
  my(t1=(a*tau+b)/(c*tau+d), NN=max(nterms(tau),nterms(t1)));
  [(c*tau+d)*thv(v,t1,NN)/ch - thv(v,tau,NN),
   (c*thv(v,t1,NN)+thd(v,t1,NN)/(c*tau+d))/ch - thd(v,tau,NN)];}
{solve3(gm,v,f1,f2) = my(a=gm[1,1],c=gm[2,1],d=gm[2,2], p0=(a-d)/(2*c));
  my(t0=p0+1/(-c/2-I*abs(c)*f1), t1=p0+1/(-c/2-I*abs(c)*f2));
  my(r0=Rval(gm,v,t0), r1=Rval(gm,v,t1));
  matsolve([1,t0,gsh(gm,t0); 0,1,gshd(gm,t0); 1,t1,gsh(gm,t1)], [r0[1],r0[2],r1[1]]~);}

gm1 = [-7,4;-16,9];             /* fold cusp 1/2, width 4 */
GG = Catalan; Z2 = Pi^2/6; PI2 = Pi^2;
NM = ["E  ","V2E","V4E","T  ","V2T","V4T"];
VV = [ainn, shft(ainn,2), shft(ainn,4), aout, shft(aout,2), shft(aout,4)];
print("=== fold cusp = the cusp 1/2 of Gamma_0(16) (x=-1/4), gamma_1=", gm1, ", width 4 ===");
print("R(tau) = alp + bet*tau + a0*g(tau);  fold-regular <=> a0=0 AND 2*alp+bet=0;  then xi=-bet/16\n");
ALP = vector(6); BET = vector(6); A0 = vector(6);
{for(j=1,6, my(z=solve3(gm1,VV[j],0.5,0.75));
  ALP[j]=z[1]; BET[j]=z[2]; A0[j]=z[3];
  print(NM[j], " a0 = ", bestappr(real(z[3]),10^20), "   (imag ", abs(imag(z[3])), ")");
  if(j<=3,
    print("     alp/G   = ", bestappr(real(z[1])/GG,10^20), "   alp/(i pi^2) = ", bestappr(imag(z[1])/PI2,10^20));
    print("     bet/G   = ", bestappr(real(z[2])/GG,10^20), "   bet/(i pi^2) = ", bestappr(imag(z[2])/PI2,10^20)),
    print("     alp/z2  = ", bestappr(real(z[1])/Z2,10^20), "   (imag ", abs(imag(z[1])), ")");
    print("     bet/z2  = ", bestappr(real(z[2])/Z2,10^20), "   (imag ", abs(imag(z[2])), ")"));
);}
print("\n=== the two linear functionals on (c1,c2,c4) ===");
print("inner  lam = -16 alp - 8 bet, in units of i*pi^2:");
{for(j=1,3, print("   ", NM[j], ": ", bestappr(imag(-16*ALP[j]-8*BET[j])/PI2, 10^20)));}
print("outer  lam = -16 alp - 8 bet, in units of zeta(2):");
{for(j=4,6, print("   ", NM[j], ": ", bestappr(real(-16*ALP[j]-8*BET[j])/Z2, 10^20),
   "  (imag ", abs(imag(-16*ALP[j]-8*BET[j])), ")"));}
print("\nxi = -bet/16 :");
{for(j=1,3, print("   ", NM[j], " : ", bestappr(real(-BET[j]/16)/GG,10^20), "*G + ",
    bestappr(imag(-BET[j]/16)/PI2,10^20), "*i*pi^2"));}
{for(j=4,6, print("   ", NM[j], " : ", bestappr(real(-BET[j]/16)/Z2,10^20), "*zeta(2)"));}
quit
