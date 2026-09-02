/* 45_eta1b.gp -- ROUTE (D), fixed precisions. */
default(parisize,10000000000);
read("40_core.gp");
{ ratfitS(S,d,p,MM) = my(Mat0, K, sol, nu, de, g);
  Mat0 = matrix(MM, 2*d+2);
  for(i=0,d, for(r=1,MM, Mat0[r,i+1] = polcoeff(S*'x^i, r-1)));
  for(i=0,d, for(r=1,MM, if(r-1==i, Mat0[r,d+2+i] = -1)));
  K = matker(Mat0*Mod(1,p));
  if(#K==0, return(0));
  sol=K[,1]; nu=sum(i=0,d, sol[d+2+i]*'T^i); de=sum(i=0,d, sol[i+1]*'T^i); g=gcd(nu,de);
  [nu/g, de/g];
}
NX = 1500;
{ EE = vector(3); SQ = vector(3); FF = vector(3); LL = vector(3); PP = vector(3);
  for(k=1,3, FF[k]=FSER(k,NX); LL[k]=LSER(k,NX); PP[k]=PSER(k,NX); SQ[k]=sqrt(PP[k]);
             EE[k]=LL[k]/('x*SQ[k]*FF[k])); }
{ ETA1(k,p) = my(e=EE[k], ps=psival(k,p), JJ=(NX-p)\p);
  vector(JJ+1, j, (polcoeff(e,p*(j-1)+p-1) - ps*polcoeff(e,j-1))/p); }
print("=== (1) identity  eta_1 = D( C(eta*h) )  mod p,  q(x)^p = q(x^p)(1+p h) ===");
MQ = 400;
{ for(k=1,3, my(xq=Setup(k,MQ)[3], qx, NT);
   qx = subst(serreverse(xq),'q,'x);
   NT = MQ - 6;
   forprime(p=2,7, if(LEV[k]%p==0,next);
     my(qxp, qxs, h, eh, JM, Ceh, coef, Dw, e1, ok=1, nt=0);
     qxp = (qx + O('x^NT))^p;
     qxs = subst(qx + O('x^(NT\p+2)), 'x, 'x^p);
     h = (qxp/qxs - 1)/p;
     eh = (EE[k] + O('x^NT))*h;
     JM = (NT-p-4)\p;
     Ceh = sum(j=0, JM, polcoeff(eh, p*j+p-1)*'x^j) + O('x^(JM+1));
     coef = Ceh * ('x*SQ[k]*FF[k]);
     Dw = deriv(coef,'x);
     e1 = ETA1(k,p);
     for(j=0, JM-2, nt++; if(Mod(e1[j+1],p) != Mod(polcoeff(Dw,j),p), ok=0));
     print("  ",NAM[k]," p=",p,": ",if(ok,"IDENTITY HOLDS","FAILS")," (",nt," coefficients tested)"))); }
print();
print("=== (2) eta_1 mod p:  is eta_1/eta rational in x? ===");
{ for(k=1,3,
   forprime(p=3,19, if(LEV[k]%p==0,next);
     my(e1=ETA1(k,p), JJ, r, res, found=0, dmax);
     JJ = #e1;
     r = Mod(1,p)*Ser(e1,'x) / (Mod(1,p)*(EE[k]+O('x^JJ)));
     dmax = (JJ-24)\3;
     for(d=0,min(12,dmax), res=ratfitS(r,d,p,3*d+20);
       if(res!=0, print("  ",NAM[k]," p=",p,": eta_1/eta = (",lift(res[1]),")/(",lift(res[2]),")"); found=1; break));
     if(!found, print("  ",NAM[k]," p=",p,": eta_1/eta NOT rational of degree <= ",min(12,dmax))))); }
print();
print("=== (2b) is eta_1 mod p a Cartier eigenvector? ===");
{ for(k=1,3, my(row=[]);
   forprime(p=3,19, if(LEV[k]%p==0,next);
     my(e1=ETA1(k,p), JJ=#ETA1(k,p), c=0, ok=1, j0=0, nt=0);
     \\ find first index with e1 nonzero mod p
     for(j=1,JJ, if(Mod(e1[j],p)!=0, j0=j; break));
     if(j0==0, row=concat(row,[[p,"eta_1 = 0 mod p"]]),
       if(p*(j0-1)+p-1+1>JJ, row=concat(row,[[p,"range too small"]]),
         c = Mod(e1[p*(j0-1)+p],p)/Mod(e1[j0],p);
         for(j=1,JJ, if(p*(j-1)+p >JJ, break); nt++;
           if(Mod(e1[p*(j-1)+p],p) != c*Mod(e1[j],p), ok=0));
         row=concat(row,[[p,lift(c),nt,if(ok,"EIGENVECTOR","not")]])));
   );
   print("  ",NAM[k],": ",row)); }
print();
print("=== (2c) lindep of eta_1 against {eta, x*eta, F'/F, P'/P, xF/l, 1} mod p ===");
{ for(k=1,3,
   forprime(p=3,17, if(LEV[k]%p==0,next);
     my(e1=ETA1(k,p), JJ, basis, Mat0, K);
     JJ = min(#e1, 200);
     basis = [EE[k], 'x*EE[k], deriv(FF[k],'x)/FF[k], deriv(PP[k],'x)/PP[k], 'x*FF[k]/LL[k], 1+O('x^JJ)];
     Mat0 = matrix(JJ, 7);
     for(r=1,JJ, Mat0[r,1] = Mod(e1[r],p));
     for(i=1,6, for(r=1,JJ, Mat0[r,i+1] = Mod(polcoeff(basis[i], r-1),p)));
     K = matker(Mat0);
     print("  ",NAM[k]," p=",p,": kernel dim = ",#K, if(#K>0, concat("  rel: ",Str(lift(K[,1]))),"")))); }
quit;
