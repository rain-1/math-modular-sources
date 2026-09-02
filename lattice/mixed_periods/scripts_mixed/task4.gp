default(parisize,"4G");
default(realprecision,900);
N=400;
intg(s)=my(r=O(x^(N+2)));for(k=0,N, r+=polcoeff(s,k)*x^(k+1)/(k+1));r;
lg=log(1-x+O(x^(N+2)));
{for(mi=1,4, my(m=[2,3,1,5][mi], D=4*m-1, a=sqrt(D+0.), th=2*atan(1/a),
  cB=th/a, cD=(th*log(D/(m+0.))-2*imag(polylog(2,exp(I*(Pi-th)))))/a,
  HA=1/sqrt(1-4*m*x+O(x^(N+2))), HB=HA*intg(HA/(1-x)), HD=HA*intg(HA*lg/(1-x)));
 print("### m=",m,"  D'=",D,"  4m=",4*m);
 printf("   c_B = %.30f   c_D = %.30f\n",cB,cD);
 print("   n | |a_n(HA)|^(1/n) | |a_n(HB)|^(1/n) | |a_n(HB)-cB a_n(HA)|^(1/n) | |a_n(HD)|^(1/n) | |a_n(HD)-cD a_n(HA)|^(1/n)");
 for(k=1,6, my(n=[25,50,100,200,300,400][k], hA=polcoeff(HA,n), hB=polcoeff(HB,n), hD=polcoeff(HD,n),
    rB=abs(hB-cB*hA), rD=abs(hD-cD*hA));
   printf("   %4d | %.14f | %.14f | %.14f | %.14f | %.14f\n",n,(hA*1.)^(1/n),abs(hB*1.)^(1/n),rB^(1/n),abs(hD*1.)^(1/n),rD^(1/n)));
);}
quit;
