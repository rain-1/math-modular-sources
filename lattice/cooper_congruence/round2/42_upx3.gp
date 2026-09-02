/* 42_upx3.gp -- route (A) continued: is U_p(x^b) mod p a rational function of x ? */
default(parisize,8000000000);
read("40_core.gp");
{ ratfit(S,t,d,MM,p) = my(cols=2*d+2, Mat0, tp, v, K, sol, nu, de);
  Mat0 = matrix(MM, cols);
  tp = 1 + O('q^(MM+3));
  v = vector(d+1);
  for(i=0,d, v[i+1] = tp; tp = tp*t);
  for(i=0,d, for(r=1,MM, Mat0[r,i+1] = polcoeff(S*v[i+1], r-1)));
  for(i=0,d, for(r=1,MM, Mat0[r,d+2+i] = -polcoeff(v[i+1], r-1)));
  K = matker(Mat0*Mod(1,p));
  if(#K==0, return(0));
  sol = K[,1];
  [sum(i=0,d, sol[d+2+i]*'T^i), sum(i=0,d, sol[i+1]*'T^i)];
}
M = 500;
print("=== is U_p(x^b) mod p rational in x?  rows s7, s10 ===");
{ for(k=1,2, my(S=Setup(k,M), x=S[3], N=LEV[k]);
   print(" row ",NAM[k]);
   forprime(p=3,13, if(N%p==0,next);
     my(xp = Mod(1,p)*x, xb = 1+O('q^(M-4)));
     for(b=1,min(p-1,3),
       xb = xb*x;
       my(nc = (M-20)\p, cf, Ub, res, found=0);
       cf = vector(nc, m, polcoeff(xb, p*m));
       Ub = Mod(1,p)*Ser(concat([0],cf),'q);
       for(d=1,16, res = ratfit(Ub,xp,d,3*d+16,p);
         if(res!=0, print("   p=",p," b=",b,": U_p(x^b) = N/D, deg<=",d,"  N=",lift(res[1]),"  D=",lift(res[2])); found=1; break));
       if(!found, print("   p=",p," b=",b,": NO rational fit in x, deg <= 16")))));
}
quit;
