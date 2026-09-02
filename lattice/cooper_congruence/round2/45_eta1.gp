/* 45_eta1.gp -- ROUTE (D).  eta_1 defined by  C(eta) = psi(p) eta + p eta_1.
   (1) verify the round-1 identity  eta_1 = D( C(eta*h) ) mod p, where q(x)^p = q(x^p)(1+p h)
       and, for w = g dx,  D(w) := d( g * x*sqrt(P)*F )  (the coefficient of dq/q is g*x*y*F).
   (2) look for an independent description of eta_1: is eta_1/eta rational in x mod p?
       is eta_1 a Cartier eigenvector?  lindep against eta, x*eta, 1/x, F'/F, P'/P, l'/l.   */
default(parisize,8000000000);
read("40_core.gp");
{ ratfitS(S,d,p) = my(MM=3*d+20, Mat0, K, sol, nu, de, g);
  Mat0 = matrix(MM, 2*d+2);
  for(i=0,d, for(r=1,MM, Mat0[r,i+1] = polcoeff(S*'x^i, r-1)));
  for(i=0,d, for(r=1,MM, if(r-1==i, Mat0[r,d+2+i] = -1)));
  K = matker(Mat0*Mod(1,p));
  if(#K==0, return(0));
  sol=K[,1]; nu=sum(i=0,d, sol[d+2+i]*'T^i); de=sum(i=0,d, sol[i+1]*'T^i); g=gcd(nu,de);
  [nu/g, de/g];
}
NX = 400;
MQ = 220;
print("=== (1) the identity  eta_1 = D(C(eta*h))  mod p ===");
{ for(k=1,3, my(R=ROWS[k],B=R[3],C=R[4], A, F, l, P, sq, e, xq, qx, row=[]);
   F=FSER(k,NX); l=LSER(k,NX); P=PSER(k,NX); sq=sqrt(P); e=l/('x*sq*F);
   xq = Setup(k,MQ)[3];
   qx = serreverse(xq);                       \\ q as a series in x  (variable q -> renamed)
   qx = subst(qx,'q,'x);
   forprime(p=3,7, if(LEV[k]%p==0,next);
     my(NT=min(NX-4, MQ-4), qxp, qxpp, h, eh, Ceh, coef, Dw, eta1, ok=1, nt=0);
     qxp  = qx + O('x^NT);
     qxpp = qxp^p;                               \\ q(x)^p
     my(qxs = subst(qxp,'x,'x^p) + O('x^NT));    \\ q(x^p)
     h = (qxpp/qxs - 1)/p;
     eh = e*h;
     Ceh = sum(j=0, (NT-p)\p, polcoeff(eh, p*j+p-1)*'x^j) + O('x^((NT-p)\p));
     coef = Ceh * ('x*sq*F);
     Dw = deriv(coef,'x);
     eta1 = sum(j=0, (NX-p)\p, ((polcoeff(e,p*j+p-1) - psival(k,p)*polcoeff(e,j))/p)*'x^j);
     for(j=0, min(20, (NT-p)\p - 3), nt++;
       if(Mod(polcoeff(eta1,j),p) != Mod(polcoeff(Dw,j),p), ok=0));
     print("  ",NAM[k]," p=",p,": ",if(ok,"IDENTITY HOLDS","FAILS"),"  (",nt," coefficients)"));
  ); }
print();
print("=== (2) is eta_1/eta mod p a rational function of x? ===");
{ for(k=1,3, my(F=FSER(k,NX), l=LSER(k,NX), P=PSER(k,NX), sq=sqrt(P), e=l/('x*sq*F));
   forprime(p=3,17, if(LEV[k]%p==0,next);
     my(eta1, r, res, found=0);
     eta1 = sum(j=0,(NX-p)\p, ((polcoeff(e,p*j+p-1)-psival(k,p)*polcoeff(e,j))/p)*'x^j) + O('x^((NX-p)\p));
     r = Mod(1,p)*eta1/(Mod(1,p)*e);
     for(d=0,10, res=ratfitS(r,d,p); if(res!=0,
        print("  ",NAM[k]," p=",p,":  eta_1/eta = ",lift(res[1]),"  /  ",lift(res[2])); found=1; break));
     if(!found, print("  ",NAM[k]," p=",p,":  eta_1/eta NOT rational of degree <= 10"))));
}
print();
print("=== (2b) is eta_1 a Cartier eigenvector?  C(eta_1) vs c*eta_1 ===");
{ for(k=1,3, my(F=FSER(k,NX), l=LSER(k,NX), P=PSER(k,NX), sq=sqrt(P), e=l/('x*sq*F), row=[]);
   forprime(p=3,17, if(LEV[k]%p==0,next);
     my(eta1, JJ, c, ok);
     JJ = (NX-p)\p;
     eta1 = vector(JJ+1, j, (polcoeff(e,p*(j-1)+p-1)-psival(k,p)*polcoeff(e,j-1))/p);
     if(eta1[1]==0, row=concat(row,[[p,"eta_1(0)=0"]]), 
       c = Mod(eta1[p],p)/Mod(eta1[1],p);      \\ candidate eigenvalue from j=0 slot
       ok=1;
       for(j=1,(JJ-p+1)\p, if(p*j+p-1+1>JJ+1,break);
         if(Mod(eta1[p*j+p],p) != c*Mod(eta1[j+1],p), ok=0));
       row=concat(row,[[p,lift(c),if(ok,"EIGEN","not eigen")]])));
   print("  ",NAM[k],": ",row)); }
print();
print("=== (2c) lindep of eta_1 against {eta, x*eta, 1/x, F'/F, P'/P, F/l} mod p ===");
{ for(k=1,3, my(F=FSER(k,NX), l=LSER(k,NX), P=PSER(k,NX), sq=sqrt(P), e=l/('x*sq*F));
   forprime(p=3,13, if(LEV[k]%p==0,next);
     my(JJ=(NX-p)\p, eta1, basis, Mat0, K);
     eta1 = vector(JJ, j, Mod((polcoeff(e,p*(j-1)+p-1)-psival(k,p)*polcoeff(e,j-1))/p, p));
     basis = [e, 'x*e, deriv(F,'x)/F, deriv(P,'x)/P, F/l*'x];
     Mat0 = matrix(JJ, 6);
     for(r=1,JJ, Mat0[r,1] = eta1[r]);
     for(i=1,5, for(r=1,JJ, Mat0[r,i+1] = Mod(polcoeff(basis[i], r-1),p)));
     K = matker(Mat0);
     print("  ",NAM[k]," p=",p,": kernel dim = ",#K, if(#K>0, concat("   relation: ",Str(lift(K[,1]))), "")));
  ); }
quit;
