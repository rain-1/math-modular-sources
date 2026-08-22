default(parisizemax, 6000000000);
default(realprecision, 400);
NV = 200; NK = 120;
LCT = vector(NK+2); LCT[1]=1; for(n=1,NK+1, LCT[n+1]=lcm(LCT[n],n));
isok(z) = my(p=lift(z)); (denominator(polcoeff(p,0))==1)&&(denominator(polcoeff(p,1))==1);
isQ(z)  = my(p=lift(z)); polcoeff(p,1)==0;
cf(z)   = my(p=lift(z)); [polcoeff(p,0),polcoeff(p,1)];
emb(z,tt,ss,sg) = my(p=lift(z)); polcoeff(p,0)+polcoeff(p,1)*(tt+sg*sqrt(tt^2+4*ss))/2;
lines = readstr(INFILE);
print("#rows read: ", #lines);
seen = Map(); nclass=0; nint=0; nnondeg=0;
{
for(li=1,#lines,
  my(f=strsplit(lines[li]," "));
  if(#f<8, next);
  my(sh=if(f[1]=="R3",3,2), D=eval(f[2]),
     a0=eval(f[3]),a1=eval(f[4]),b0=eval(f[5]),b1=eval(f[6]),c0=eval(f[7]),c1=eval(f[8]),
     tt,ss,k1,k2,k3,k4,kk);
  if(D%4==1, tt=1; ss=(D-1)/4, tt=0; ss=D);
  k1=[sh,D,a0,a1,b0,b1,c0,c1];
  k2=[sh,D,-a0,-a1,-b0,-b1,c0,c1];
  k3=[sh,D,a0+a1*tt,-a1,b0+b1*tt,-b1,c0+c1*tt,-c1];
  k4=[sh,D,-(a0+a1*tt),a1,-(b0+b1*tt),b1,c0+c1*tt,-c1];
  kk=vecsort([k1,k2,k3,k4])[1];
  if(mapisdefined(seen,kk), next);
  mapput(seen,kk,1); nclass++;
  my(ww=Mod(y,y^2-tt*y-ss), AA=a0+a1*ww, BB=b0+b1*ww, CC=c0+c1*ww);
  my(deg = if(sh==3, AA^2==CC, AA^2==4*CC) || CC==0 || (AA==0 && BB==0));
  my(rat = (a1==0)&&(b1==0)&&(c1==0));
  my(twq = (AA!=0) && isQ(BB/AA) && isQ(CC/AA^2));
  my(uu=vector(NV+2), r=sh, bad=-1);
  uu[1]=1; uu[2]=BB;
  if(sh==3,
    for(n=1,NV, uu[n+2]=((2*n+1)*(AA*n^2+AA*n+BB)*uu[n+1]-CC*n^3*uu[n])/(n+1)^3),
    for(n=1,NV, uu[n+2]=((AA*n^2+AA*n+BB)*uu[n+1]-CC*n^2*uu[n])/(n+1)^2));
  for(n=0,NV, if(!isok(uu[n+1]), bad=n; break));
  if(bad>=0, print("FALSEPOS ",lines[li]," first-nonintegral n=",bad); next);
  nint++; if(!deg, nnondeg++);
  my(cb=vector(NK+2)); cb[1]=0; cb[2]=1;
  if(sh==3,
    for(n=1,NK, cb[n+2]=((2*n+1)*(AA*n^2+AA*n+BB)*cb[n+1]-CC*n^3*cb[n])/(n+1)^3),
    for(n=1,NK, cb[n+2]=((AA*n^2+AA*n+BB)*cb[n+1]-CC*n^2*cb[n])/(n+1)^2));
  my(sk=-1, skf=-1);
  for(kx=0,8, my(fl=-1); for(n=1,NK, if(!isok(LCT[n+1]^kx*cb[n+1]), fl=n; break));
              if(fl<0, sk=kx; break));
  if(sk>0, my(fl=-1); for(n=1,NK, if(!isok(LCT[n+1]^(sk-1)*cb[n+1]), fl=n; break)); skf=fl);
  my(l1=vector(2), l2=vector(2), cpx=0);
  for(i=1,2, my(sg=if(i==1,1,-1), aa=emb(AA,tt,ss,sg), cc=emb(CC,tt,ss,sg), ac,dd,r1,r2);
     if(sh==3, ac=aa; dd=aa^2-cc, ac=aa/2; dd=aa^2/4-cc);
     if(dd>=0, r1=abs(ac+sqrt(dd)); r2=abs(ac-sqrt(dd)), r1=sqrt(abs(cc)); r2=r1; cpx=1);
     if(r1<r2, my(z=r1); r1=r2; r2=z);
     l1[i]=r1; l2[i]=r2);
  my(ku=if(sk>0,sk,sh));
  my(mA=-(log(l2[1])+log(l1[2]))/2-ku, mB=-(log(l2[2])+log(l1[1]))/2-ku);
  my(mI=max(mA,mB), v0=if(mA>=mB,1,2), mII=-(log(l2[1])+log(l2[2]))/2-ku);
  print("ROW|",sh,"|",D,"|",a0,"|",a1,"|",b0,"|",b1,"|",c0,"|",c1,"|",sk,"|",skf,"|",
        deg,"|",rat,"|",twq,"|",precision(l1[1],20),"|",precision(l2[1],20),"|",precision(l1[2],20),"|",precision(l2[2],20),"|",precision(mI,20),"|",v0,"|",precision(mII,20),"|",cpx);
);
}
print("SUMMARY classes=",nclass," integral=",nint," nondegenerate=",nnondeg);
quit;
