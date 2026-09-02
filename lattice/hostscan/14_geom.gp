/* exact-recognised special-point images for the 13 Fricke families */
default(realprecision, 60);
uev(dv, r, t) = prod(j=1, #dv, eta(dv[j]*t, 1)^r[j]);
cusplist(n) = my(res=List()); fordiv(n, c, my(g=gcd(c,n/c), seen=List()); for(a=1, max(g,1), if(gcd(a,c)==1 || c==1, if(gcd(a,g)==1, my(new=1); for(k=1,#seen, if((a-seen[k])%g==0, new=0; break)); if(new, listput(seen,a); listput(res,[a,c])))))); Vec(res);
cuspmat(a,c) = my(g=bezout(a,c)); [a, -g[2]; c, g[1]];
ell2(n) = my(res=List()); for(d=0, n-1, if((d^2+1)%n==0, listput(res, (-d+I)/n))); Vec(res);
ell3(n) = my(res=List()); for(d=0, n-1, if((d^2-d+1)%n==0, listput(res, ((1-2*d)+I*sqrt(3))/(2*n)))); Vec(res);
rec(v) = my(re=bestappr(real(v),10^12), im=bestappr(imag(v),10^12)); if(abs(imag(v))<1e-30, Str(re), Str(re," + (",im,")i"));
Y = 90.0;
{PARS = [[5,[1,5],[-6,6],125], [6,[1,2,3,6],[-5,1,-1,5],72], [7,[1,7],[-4,4],49], [8,[1,2,4,8],[-4,2,-2,4],32], [9,[1,3,9],[-3,0,3],27], [10,[1,2,5,10],[-3,1,-1,3],20], [12,[1,2,3,4,6,12],[-3,2,1,-1,-2,3],12], [13,[1,13],[-2,2],13], [16,[1,2,4,8,16],[-2,1,0,-1,2],8], [18,[1,2,3,6,9,18],[-2,1,1,-1,-1,2],6]];}
{for(h=1,#PARS, my(nn=PARS[h][1], dv=PARS[h][2], r=PARS[h][3], C=PARS[h][4]);
  print("N=",nn," C=",C," u=eta^",r);
  my(cl = cusplist(nn));
  for(i=1,#cl, my(a=cl[i][1], c=cl[i][2], m=cuspmat(a,c));
    my(t=(m[1,1]*(I*Y)+m[1,2])/(m[2,1]*(I*Y)+m[2,2])); my(v=uev(dv,r,t));
    if(abs(v)<1e-30, print("   cusp ",a,"/",c,"  u = 0      (x=0)"),
      if(abs(v)>1e30, print("   cusp ",a,"/",c,"  u = oo     (x=0, W_N image)"),
        print("   cusp ",a,"/",c,"  u = ",rec(v),"   1/u = ",rec(1/v)))));
  my(e2=ell2(nn)); for(i=1,#e2, my(v=uev(dv,r,e2[i])); print("   ell2       u = ",rec(v),"   1/u = ",rec(1/v)));
  my(e3=ell3(nn)); for(i=1,#e3, my(v=uev(dv,r,e3[i])); print("   ell3       u = ",rec(v),"   1/u = ",rec(1/v)));
  print("   W_N fixed points (the FOLDS): u = +-1/sqrt(",C,") = +-",1/sqrt(C));
);}
quit;
