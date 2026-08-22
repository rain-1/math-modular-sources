\\ Companion b_n of a root row whose fitted recurrence has degree 2 in n
\\ (i.e. the differential operator L' = sum_j u^j q_j(theta) is SECOND ORDER).
\\ Solve L'(B) = c*u with b_0=0, b_1=1 and test d_n^k b_n in Z  (Theorem R3 predicts k=2).
read("lattice/root_rows/lib.gp");
default(realprecision,40);
readrow(name)={readvec(Str("lattice/root_rows/rows_",name,".txt"))};

companion(name,w,o,NB)={
 my(A=readrow(name), rr=rootrow(A,w,#A-1), a=rr[2], N=#a-1, lam=rr[1]);
 my(r=fitrec(a,o,2));
 print("\n#### ",name,"  w=",w,"  lambda=",lam,"  recurrence order ",o,", degree 2");
 my(qs=vector(o+1,j,subst(r[j],'n,'n+(j-1))));
 print("  q_0(th) = ",qs[1],"   (differential operator is order ",poldegree(qs[1]),")");
 \\ extend a_n to NB using the recurrence
 my(av=vector(NB+1)); for(i=0,min(N,NB), av[i+1]=a[i+1]);
 for(m=N+1,NB, my(c=subst(r[1],'n,m)); if(c==0,print("  leading vanishes at ",m);break);
   av[m+1]=-sum(j=1,o, subst(r[j+1],'n,m)*av[m-j+1])/c);
 my(bv=vector(NB+1)); bv[1]=0; bv[2]=1;
 for(m=2,NB, my(c=subst(qs[1],'n,m)); if(c==0,print("  q_0 vanishes at ",m);break);
   bv[m+1] = -sum(j=1,o, if(m-j>=0, subst(qs[j+1],'n,m-j)*bv[m-j+1], 0))/c);
 print("  a_0..a_4 = ",vector(5,j,av[j]));
 print("  b_1..b_5 = ",vector(5,j,bv[j+1]));
 my(dn=1,K=[1,1,1,1,1]);
 for(n=1,NB, dn=lcm(dn,n); for(kk=1,5, if(denominator(dn^kk*bv[n+1])!=1, K[kk]=0)));
 print("  d_n^k b_n in Z (k=1..5), n<=",NB,": ",K);
 my(rt=charroots(r), md=vecsort(vector(#rt,i,abs(rt[i]))));
 print("  |char roots| = ",md);
 print("  lambda_1 = ",md[#md],"   next distinct = ",md[#md-1]);
 my(xi=bv[NB+1]/av[NB+1]*1.0);
 print("  b_n/a_n at n=",NB," = ",xi,"   diff to n-1: ",abs(xi-bv[NB]/av[NB])*1.0);
 [lam,r,md,xi]};
C1=companion("zeta5_L16",4,14,400);
C2=companion("zeta7_L24",6,10,400);
\q
