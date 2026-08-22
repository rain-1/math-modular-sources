default(parisizemax, 6000000000);
default(realprecision, 60);
isok(z)=my(p=lift(z)); (denominator(polcoeff(p,0))==1)&&(denominator(polcoeff(p,1))==1);
cf(z)=my(p=lift(z)); [polcoeff(p,0),polcoeff(p,1)];
LCT=vector(202); LCT[1]=1; for(n=1,201, LCT[n+1]=lcm(LCT[n],n));
\\ rows: [D, a0,a1, b0,b1, c0,c1]
rows = [[5,-110,-179,-30,-49,1,1],[5,-13,-28,-3,-8,13,-8],[5,-179,110,-49,30,1,0],[5,-28,-41,-8,-11,5,-3],[5,41,69,11,19,2,-1]];
{
for(i=1,#rows,
  my(v=rows[i], D=v[1], a0=v[2],a1=v[3],b0=v[4],b1=v[5],c0=v[6],c1=v[7], tt,ss);
  if(D%4==1, tt=1; ss=(D-1)/4, tt=0; ss=D);
  my(ww=Mod(y,y^2-tt*y-ss), AA=a0+a1*ww, BB=b0+b1*ww, CC=c0+c1*ww);
  my(sq=2*ww-tt);  \\ sqrt(D) when D=1 mod4 (w=(1+sqrtD)/2 -> 2w-1); when tt=0, 2w = 2 sqrtD -- fix below
  print("=========================================================");
  print("D=",D,"  A = ",a0," + ",a1,"*w   B = ",b0," + ",b1,"*w   C = ",c0," + ",c1,"*w   (w=(1+sqrt",D,")/2)");
  print("  in a+b*sqrt",D," form:  A = ",[a0+a1/2,a1/2],"  B = ",[b0+b1/2,b1/2],"  C = ",[c0+c1/2,c1/2]);
  print("  N(A)=",norm(AA)," N(B)=",norm(BB)," N(C)=",norm(CC));
  my(u=vector(202)); u[1]=1; u[2]=BB;
  for(n=1,200, u[n+2]=((2*n+1)*(AA*n^2+AA*n+BB)*u[n+1]-CC*n^3*u[n])/(n+1)^3);
  my(bad=-1); for(n=0,200, if(!isok(u[n+1]), bad=n; break));
  print("  integral for 0<=n<=200 : ", if(bad<0,"YES",Str("NO, fails at n=",bad)));
  print("  first ten u_n  (as [x,y] meaning x + y*sqrt",D,"):");
  for(j=0,9, my(z=cf(u[j+1])); print("    u_",j," = [",z[1]+z[2]/2,", ",z[2]/2,"]"));
  my(cb=vector(122)); cb[1]=0; cb[2]=1;
  for(n=1,120, cb[n+2]=((2*n+1)*(AA*n^2+AA*n+BB)*cb[n+1]-CC*n^3*cb[n])/(n+1)^3);
  for(kx=1,4, my(fl=-1); for(n=1,120, if(!isok(LCT[n+1]^kx*cb[n+1]), fl=n; break));
      print("    companion: lcm^",kx," integral to n=120? ", if(fl<0,"YES",Str("NO, first fail n=",fl))));
  my(l1=vector(2), l2=vector(2), rt=vector(2));
  for(q=1,2, my(sg=if(q==1,1,-1), aa=(a0+a1*(tt+sg*sqrt(tt^2+4*ss))/2), cc=(c0+c1*(tt+sg*sqrt(tt^2+4*ss))/2), dd, r1,r2);
     dd=aa^2-cc;
     if(dd>=0, r1=aa+sqrt(dd); r2=aa-sqrt(dd), r1=0; r2=0);
     rt[q]=[r1,r2];
     r1=abs(r1); r2=abs(r2); if(r1<r2, my(z=r1);r1=r2;r2=z); l1[q]=r1; l2[q]=r2);
  print("  place v1 (sqrt",D,"->+): roots ",rt[1],"   |l1|=",l1[1]," |l2|=",l2[1]);
  print("  place v2 (sqrt",D,"->-): roots ",rt[2],"   |l1|=",l1[2]," |l2|=",l2[2]);
  my(m1=-(log(l2[1])+log(l1[2]))/2-3, m2=-(log(l2[2])+log(l1[1]))/2-3);
  print("  ModeI (v0=v1) = ",m1,"   ModeI (v0=v2) = ",m2);
  print("  ModeI margin  = ",max(m1,m2)," (v0=v",if(m1>=m2,1,2),")   geom mean of radii = ",exp(max(m1,m2)+3));
  print("  N(lambda2) = ",l2[1]*l2[2],"   sqrt = ",sqrt(l2[1]*l2[2]),"   ModeII margin = ",-(log(l2[1])+log(l2[2]))/2-3);
);
}
quit;
