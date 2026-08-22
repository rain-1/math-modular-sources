/* 17_conditions.gp -- from the exact (alp,bet) at each cusp, print the linear
   fold-regularity condition on the inner coordinates (c1,c2,c4); only the fold
   cusp 1/2 can reproduce the condition 8c2+c4=0 verified in CATALAN_TWO_CLASSES.md */
default(parisizemax, 8*10^9);
default(realprecision, 130);
PRECDIG = 130; NTOP = 6000;
chim(n) = if(n%2==0, 0, if(n%4==1, 1, -1));
ainn = vector(NTOP, n, sumdiv(n, d, chim(n/d)*d^2));
shft(v,dd) = vector(#v, n, if(n%dd==0, v[n/dd], 0));
thv(v,tau,NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n^2*p)); s;
thd(v,tau,NN) = my(qq=exp(2*Pi*I*tau), s=0., p=1.); for(n=1,NN, p*=qq; if(v[n]!=0, s += v[n]/n*p)); 2*Pi*I*s;
nterms(tau) = my(y=imag(tau)); ceil((PRECDIG+25)*log(10)/(2*Pi*y));
{ab(gm, v) = my(a=gm[1,1],b=gm[1,2],c=gm[2,1],d=gm[2,2],ch=chim(d), p0=(a-d)/(2*c));
  my(t0 = p0 + 1/(-c/2 - I*abs(c)/2), t1=(a*t0+b)/(c*t0+d));
  my(NN = max(nterms(t0),nterms(t1)));
  my(R = (c*t0+d)*thv(v,t1,NN)/ch - thv(v,t0,NN));
  my(Rp = (c*thv(v,t1,NN) + thd(v,t1,NN)/(c*t0+d))/ch - thd(v,t0,NN));
  [R - Rp*t0, Rp];}
{CU = [["0  ",[1,0;-16,1]], ["1/2",[-7,4;-16,9]], ["1/4",[-3,1;-16,5]],
       ["3/4",[-11,9;-16,13]], ["1/8",[-7,1;-64,9]]];}
VS = [ainn, shft(ainn,2), shft(ainn,4)];
{for(k=1,#CU, my(gm=CU[k][2], c=gm[2,1], d=gm[2,2], co=vector(3));
  for(j=1,3, my(z=ab(gm,VS[j])); co[j] = c*z[1] - (d-1)*z[2]);
  my(nz=0); for(j=1,3, if(abs(co[j])>1e-100, nz=j));
  print("cusp ", CU[k][1], "  condition coefficients (c1,c2,c4) = ",
     if(nz==0, "VACUOUS", [bestappr(co[1]/co[nz],10^18), bestappr(co[2]/co[nz],10^18), bestappr(co[3]/co[nz],10^18)]));)}
quit
