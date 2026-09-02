/* where does the pole pair of x sit?  cusp values of u for the 8 deg>1 / elliptic families */
default(realprecision, 50);
uev(dv, r, t) = prod(j=1, #dv, eta(dv[j]*t, 1)^r[j]);
cusplist(n) = my(res=List()); fordiv(n, c, my(g=gcd(c,n/c), seen=List()); for(a=1, max(g,1), if(gcd(a,c)==1 || c==1, if(gcd(a,g)==1, my(new=1); for(k=1,#seen, if((a-seen[k])%g==0, new=0; break)); if(new, listput(seen,a); listput(res,[a,c])))))); Vec(res);
cuspmat(a,c) = my(g=bezout(a,c)); [a, -g[2]; c, g[1]];
rec(v) = my(re=bestappr(real(v),10^10), im=bestappr(imag(v),10^10)); if(abs(imag(v))<1e-25, Str(re), Str(re," + (",im,")i"));
Y=80.0;
doit(nn,dv,r,C,B,nm) = print("== ",nm,"   poles of x: roots of 1+",B,"u+",C,"u^2  = ",polroots(C*x^2+B*x+1)); my(cl=cusplist(nn)); for(i=1,#cl, my(a=cl[i][1], c=cl[i][2], m=cuspmat(a,c)); my(t=(m[1,1]*(I*Y)+m[1,2])/(m[2,1]*(I*Y)+m[2,2])); my(v=uev(dv,r,t)); if(abs(v)<1e-25, print("   cusp ",a,"/",c," u=0"), if(abs(v)>1e25, print("   cusp ",a,"/",c," u=oo"), print("   cusp ",a,"/",c," u=",rec(v)))));
doit(6,[1,2,3,6],[-4,-4,4,4],81,14,"N=6 C=81 B=14 (zeta(2)/8, meromorphic)");
doit(6,[1,2,3,6],[-6,6,-6,6],64,20,"N=6 C=64 B=20 (L(2,chi-3)/4, meromorphic)");
doit(8,[1,2,4,8],[-8,16,-16,8],16,24,"N=8 C=16 B=24 (G/4, meromorphic)");
doit(10,[1,2,5,10],[-2,-2,2,2],25,6,"N=10 C=25 B=6 (Cooper s10, meromorphic)");
doit(12,[1,2,3,4,6,12],[-4,4,4,-4,-4,4],9,10,"N=12 C=9 B=10 (Domb, HOLOMORPHIC)");
doit(12,[1,2,3,4,6,12],[-12,24,12,-12,-24,12],1,34,"N=12 C=1 B=34 (meromorphic)");
doit(18,[1,2,3,6,9,18],[-6,6,12,-12,-6,6],1,14,"N=18 C=1 B=14 (Cooper s18, meromorphic)");
quit;
