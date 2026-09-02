default(parisize,"4G");
default(realprecision,900);
N=400;
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
lg=log(1-x+O(x^(N+2)));
{for(mi=1,2, my(m=[2,3][mi], D=4*m-1, a=sqrt(D+0.), th=2*atan(1/a),
  cB=th/a, cD=(th*log(D/(m+0.))-2*imag(polylog(2,exp(I*(Pi-th)))))/a,
  HA=1/sqrt(1-4*m*x+O(x^(N+2))), HB=HA*intg(HA/(1-x)), HD=HA*intg(HA*lg/(1-x)));
 print("### m=",m,"  |r_n| for r=H_B-c_B H_A and H_D-c_D H_A, and n*|r_n|");
 for(k=1,5, my(n=[50,100,200,300,400][k], hA=polcoeff(HA,n),
    rB=abs(polcoeff(HB,n)-cB*hA), rD=abs(polcoeff(HD,n)-cD*hA));
   printf("   n=%4d  |rB_n|=%.6e  n*|rB_n|=%.8f   |rD_n|=%.6e  n*|rD_n|=%.8f\n",n,rB,n*rB,rD,n*rD)));}
quit;
