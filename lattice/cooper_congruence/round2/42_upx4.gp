/* 42_upx4.gp -- route (A): the explicit mod-p formula for U_p on the x-line.
   PROVED (see FINDINGS): with W = dq/q / dx = 1/(x y F), y=sqrt(P), and
        H := P^((p-1)/2) * F^(p-1) = (yF)^(p-1)   in Z[[x]],
   one has, mod p,
        U_p(x^b) = x * sum_{j>=0} H_{pj+p-b} x^j          (1 <= b <= p-1)
        U_p(1)   = 1 + sum_{j>=1} H_{pj} x^j              (so U_p(1)=1 forces H_{pj}=0, j>=1)
        Xi|U_p   = sum_{j>=1} [x^{pj}]( l(x) H(x) ) x^j
   and the target congruence becomes the CLOSED FINITE STATEMENT
        comp_0( l * H ) = psi(p) * l   (mod p),   comp_0(S) := sum_j [x^{pj}]S x^j.
   H is W_N-invariant (F|W_N=-F, p-1 even), so H is a RATIONAL FUNCTION of x mod p;
   its denominator is the supersingular locus.  All of this is checked here.        */
default(parisize,8000000000);
read("40_core.gp");
{ ratfitS(S,d,p) = my(MM=3*d+18, Mat0, K, sol);
  Mat0 = matrix(MM, 2*d+2);
  for(i=0,d, for(r=1,MM, Mat0[r,i+1] = polcoeff(S*'x^i, r-1)));
  for(i=0,d, for(r=1,MM, if(r-1==i, Mat0[r,d+2+i] = -1)));
  K = matker(Mat0*Mod(1,p));
  if(#K==0, return(0));
  sol = K[,1];
  my(nu=sum(i=0,d, sol[d+2+i]*'T^i), de=sum(i=0,d, sol[i+1]*'T^i), g=gcd(nu,de));
  [nu/g, de/g];
}
NX = 420;
print("=== H = (sqrt(P) F)^(p-1) mod p as a rational function of x ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], F=FSER(k,NX), P=1-2*B*'x+(B^2-4*C)*'x^2, sq);
   print(" row ",NAM[k]);
   forprime(p=3,17, if(LEV[k]%p==0,next);
     my(H = Mod(1,p)*(P^((p-1)/2))*(F^(p-1)), res, found=0);
     for(d=0,14, res=ratfitS(H,d,p); if(res!=0,
        print("   p=",p,":  H = N/D,  deg N=",poldegree(res[1]),"  deg D=",poldegree(res[2]),
              "   D=",lift(res[2]/polcoeff(res[2],poldegree(res[2])))); found=1; break));
     if(!found, print("   p=",p,":  H not rational of degree <= 14"))));
}
print();
print("=== check  [x^{pj}] H = 0  for j>=1  (equivalent to U_p(1)=1) ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], F=FSER(k,NX), P=1-2*B*'x+(B^2-4*C)*'x^2, row=[]);
   forprime(p=3,29, if(LEV[k]%p==0,next);
     my(H=Mod(1,p)*(P^((p-1)/2))*(F^(p-1)), bad=0);
     for(j=1,(NX-8)\p, if(polcoeff(H,p*j)!=0, bad++));
     row=concat(row,[[p,(NX-8)\p,bad]]));
   print("  ",NAM[k]," [p,#tests,#nonzero]: ",row)); }
print();
print("=== check  comp_0(l*H) = psi(p)*l  (mod p)  -- the closed form of (S) mod p ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], F=FSER(k,NX), l=LSER(k,NX), P=1-2*B*'x+(B^2-4*C)*'x^2, row=[]);
   forprime(p=3,29, if(LEV[k]%p==0,next);
     my(H=Mod(1,p)*(P^((p-1)/2))*(F^(p-1)), S=Mod(1,p)*l*H, ps=psival(k,p), bad=0, nt=0);
     for(j=1,(NX-8)\p, nt++; if(polcoeff(S,p*j) != ps*Mod(polcoeff(l,j),p), bad++));
     row=concat(row,[[p,nt,bad]]));
   print("  ",NAM[k]," [p,#tests,#fail]: ",row)); }
print();
print("=== cross-check formula U_p(x^b) = x*sum_j H_{pj+p-b}x^j against q-expansions (s7) ===");
{ my(MQ=300, S=Setup(1,MQ), xq=S[3], R=ROWS[1],B=R[3],C=R[4], F=FSER(1,NX), P=1-2*B*'x+(B^2-4*C)*'x^2);
  forprime(p=3,13, if(7%p==0,next);
    my(H=Mod(1,p)*(P^((p-1)/2))*(F^(p-1)), xb=1+O('q^MQ), ok=1);
    for(b=1,min(p-1,3), xb=xb*xq;
      my(pred = Mod(1,p)*'x*sum(j=0,(NX-10)\p, polcoeff(H,p*j+p-b)*'x^j) + O('x^60), predq, nc);
      predq = subst(pred, 'x, Mod(1,p)*xq + O('q^59));
      nc = min(50, (MQ-10)\p);
      for(m=1,nc, if(Mod(polcoeff(xb,p*m),p) != polcoeff(predq,m), ok=0)));
    print("  p=",p,": ",if(ok,"FORMULA CONFIRMED","MISMATCH")));
}
quit;
