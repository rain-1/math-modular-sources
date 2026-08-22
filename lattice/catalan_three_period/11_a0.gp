/* Solve for the constant term a_0 of the outer Eisenstein series from the
   requirement that R(tau) - a_0*g(tau) be affine in tau,
   g(tau) = 2 pi^2 [ chi(d)^{-1}(c tau+d)(g tau)^2 - tau^2 ].              */
default(parisizemax, 8*10^9);
default(realprecision, 130);
PRECDIG = 130;
NTOP = 6000;
chim(n) = if(n%2==0, 0, if(n%4==1, 1, -1));
ainn = vector(NTOP, n, sumdiv(n, d, chim(n/d)*d^2));
aout = vector(NTOP, n, sumdiv(n, d, chim(d)*d^2));
shft(v, dd) = vector(#v, n, if(n%dd==0, v[n/dd], 0));
thv(v, tau, NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n^2*p)); s;
thd(v, tau, NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n*p));  2*Pi*I*s;
nterms(tau) = my(y=imag(tau)); ceil((PRECDIG+25)*log(10)/(2*Pi*y));
gsh(gm, tau) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2], ch=chim(d)); 2*Pi^2*((a*tau+b)^2/(ch*(c*tau+d)) - tau^2);
gshd(gm, tau)= my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2], ch=chim(d)); 2*Pi^2*((2*a*(a*tau+b)*(c*tau+d) - c*(a*tau+b)^2)/(ch*(c*tau+d)^2) - 2*tau);
{Rval(gm, v, tau) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2], ch=chim(d));
  my(t1=(a*tau+b)/(c*tau+d), NN=max(nterms(tau),nterms(t1)));
  [(c*tau+d)*thv(v,t1,NN)/ch - thv(v,tau,NN),
   (c*thv(v,t1,NN) + thd(v,t1,NN)/(c*tau+d))/ch - thd(v,tau,NN)];}
{solve3(gm, v, f1, f2) = my(a=gm[1,1],c=gm[2,1],d=gm[2,2]);
  my(p0=(a-d)/(2*c));
  my(t0=p0+1/(-c/2 - I*abs(c)*f1), t1=p0+1/(-c/2 - I*abs(c)*f2));
  my(r0=Rval(gm,v,t0), r1=Rval(gm,v,t1));
  my(mat=[1,t0,gsh(gm,t0); 0,1,gshd(gm,t0); 1,t1,gsh(gm,t1)]);
  matsolve(mat, [r0[1],r0[2],r1[1]]~);}
{CU = [["0  ",[1,0;-16,1]], ["1/2",[-7,4;-16,9]], ["1/4",[-3,1;-16,5]], ["3/4",[-11,9;-16,13]], ["1/8",[-7,1;-64,9]]];}
{VS = [["T  ",aout],["V2T",shft(aout,2)],["V4T",shft(aout,4)],["E  ",ainn],["V2E",shft(ainn,2)]];}
{for(k=1,#CU, my(cu=CU[k]);
  print("\n--- cusp ", cu[1], "  gamma=", cu[2]);
  for(j=1,#VS, my(z=solve3(cu[2], VS[j][2], 0.5, 0.75));
    print("   ", VS[j][1], "  alp=", z[1], "  bet=", z[2], "  a0=", z[3])));}
quit
