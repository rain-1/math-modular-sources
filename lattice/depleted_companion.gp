default(parisizemax,2000000000);
MM=420; default(seriesprecision,MM+4);
et(k)=eta(q^k+O(q^(MM+4)));
tC = q*et(1)^4*et(6)^8/(et(2)^8*et(3)^4);
FC = et(2)^6*et(3)/(et(1)^3*et(6)^2);
tt=serreverse(tC); Ft=subst(FC,q,tt); AA=vector(MM-2,n,polcoeff(Ft,n));
src = FC*(q*deriv(tC,q)); cs=vector(MM,m,polcoeff(src,m)); print("src: ",cs[1..8]);
depc(pp)={my(Th=sum(m=1,MM,if(pp>1&&m%pp==0,0,cs[m]/m^2)*q^m)); my(B=subst(FC*Th,q,tt)); vector(MM-2,n,polcoeff(B,n))};
fl=depc(1);
{for(i=1,3,my(pp=[5,7,11][i],b=depc(pp));
 print("p=",pp,"  depleted: ",vector(4,k,my(n=100*k);valuation(b[n]/AA[n]-b[n-1]/AA[n-1],pp)),"   full: ",vector(4,k,my(n=100*k);valuation(fl[n]/AA[n]-fl[n-1]/AA[n-1],pp))));}
\q
