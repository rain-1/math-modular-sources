\\ Shared library for the root-row census.

fitrec(A,ord,d)={my(M=#A,rows=[]);
 for(n=ord+1,M-1, my(r=[]); for(j=0,ord, for(e=0,d, r=concat(r,[n^e*A[n+1-j]]))); rows=concat(rows,[r]));
 my(Mx=matrix(#rows,(ord+1)*(d+1),i,j,rows[i][j]),K=matker(Mx));
 if(#K!=1,return(0)); my(v=K[,1]); v=v/content(v);
 vector(ord+1,j,Polrev(vector(d+1,e,v[(j-1)*(d+1)+e]),'n))};

\\ minimal order/degree recurrence search
minrec(A,maxord=6,maxdeg=6)={
 for(o=1,maxord, for(d=0,maxdeg,
   if((o+1)*(d+1)+8 <= #A-o,
     my(r=fitrec(A,o,d)); if(type(r)=="t_VEC",
       \\ verify on all available terms
       my(ok=1); for(n=o+1,#A-1, my(s=sum(j=0,o, subst(r[j+1],'n,n)*A[n+1-j])); if(s!=0,ok=0;break));
       if(ok, return([o,d,r]))))));
 0};

\\ w-th root with minimal lambda; returns [lambda, vector a_n, e_p data]
rootrow(A,w,N)={my(F=sum(k=0,N,A[k+1]*'x^k)+O('x^(N+1)), g, lam=1, ep=[]);
 g = F^(1/w);
 my(gc=vector(N+1,i,polcoeff(g,i-1)));
 my(pf=factor(w)[,1]);
 for(i=1,#pf, my(p=pf[i], m=0);
   for(n=1,N, my(v=valuation(gc[n+1],p)); if(gc[n+1]!=0 && v< -m*n, m=ceil(-v/n)));
   \\ minimal exponent m with p^(m n) gc[n] integral at p
   my(mm=0); for(n=1,N, if(gc[n+1]!=0, my(v=valuation(gc[n+1],p)); if(v<0, mm=max(mm,ceil(-v/n)))));
   lam *= p^mm; ep=concat(ep,[[p,mm]]));
 my(a=vector(N+1,i,lam^(i-1)*gc[i]));
 [lam,a,ep]};

\\ e_p = min_{n>=1} v_p(A_n)
emin(A,p,N)={my(e=oo); for(n=1,N, if(A[n+1]!=0, e=min(e,valuation(A[n+1],p)))); e};

\\ characteristic roots from a 3-term recurrence  p2(n)a_{n+1}+p1(n)a_n+p0(n)a_{n-1}=0
charroots(r)={my(o=#r-1, d=max(poldegree(r[1]),0));
 for(j=1,#r, d=max(d,poldegree(r[j])));
 my(c=vector(#r,j, polcoeff(r[j],d,'n)));
 \\ a_{n+1-j} ~ lam^{-j} : sum_j c_j lam^{o-j} = 0 in lam
 my(P=sum(j=1,#r, c[j]*'L^(o+1-j)));
 polroots(P)};

dlcm(n)={my(d=1);for(k=1,n,d=lcm(d,k));d};
