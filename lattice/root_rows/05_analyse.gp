read("lattice/root_rows/lib.gp");
default(realprecision, 40);

readrow(name)={my(v=readvec(Str("lattice/root_rows/rows_",name,".txt"))); v};

analyse(name, w, NBIG)={
 my(A=readrow(name), N=#A-1);
 print("\n############ ",name,"   w=",w,"   (",N+1," terms) ############");
 print("A_0..A_6 = ",vector(7,j,A[j]));
 my(intA=1); for(n=0,N, if(denominator(A[n+1])!=1, intA=0;break));
 print("A_n integral: ",intA);
 my(pf=factor(w)[,1]);
 for(i=1,#pf, print("  e_",pf[i]," = min v_p(A_n) = ",emin(A,pf[i],N),
    "   graded bound lambda_p = ",pf[i]^max(0,valuation(w,pf[i])+1-emin(A,pf[i],N))));
 my(rr=rootrow(A,w,N), lam=rr[1], a=rr[2]);
 print("  MINIMAL lambda (measured, n<=",N,") = ",lam,"     [Thm R1 bound w*rad(w) = ",w*factorback(pf),"]");
 print("  a_0..a_6 = ",vector(7,j,a[j]));
 my(ok=1); for(n=0,N, if(denominator(a[n+1])!=1, ok=0; print("  a_n NOT integral at n=",n); break));
 print("  a_n integral to n=",N,": ",ok);
 my(mr=minrec(a,5,6));
 if(type(mr)=="t_INT", print("  *** NO recurrence of order<=5, degree<=6 ***"); return(0));
 print("  minimal recurrence order=",mr[1],"  degree=",mr[2]);
 print("     ",mr[3]);
 my(rts=charroots(mr[3]));
 print("  characteristic roots: ",rts);
 if(mr[1]!=2, print("  *** NOT second order -- root row is not a rank-2 object in this coordinate ***"); return(0));
 \\ extend a_n and build b_n
 my(r=mr[3], av=vector(NBIG+1), bv=vector(NBIG+1));
 for(i=0,N, av[i+1]=a[i+1]);
 bv[1]=0; bv[2]=1;
 my(p2=r[1],p1=r[2],p0=r[3]);
 for(n=1,NBIG-1,
   my(c2=subst(p2,'n,n+1), c1=subst(p1,'n,n+1), c0=subst(p0,'n,n+1));
   \\ recurrence indexing: sum_j r[j+1](n) a_{n+1-j} = 0  =>  r1(n)a_n + r2(n)a_{n-1}+r3(n)a_{n-2}=0
   \\ so with m=n+1:  r1(m)a_m + r2(m)a_{m-1} + r3(m)a_{m-2}=0 ; solve for a_m needs r1(m)!=0
   0);
 \\ rewrite: fitrec gives sum_{j=0}^{o} r[j+1](n) * A[n+1-j] = 0 for n>=o+1, i.e.
 \\   r[1](n)a_n + r[2](n)a_{n-1} + r[3](n)a_{n-2} = 0.
 for(m=N+1,NBIG,
   my(c=subst(r[1],'n,m)); if(c==0, print("  leading coeff vanishes at n=",m); break);
   av[m+1] = -(subst(r[2],'n,m)*av[m] + subst(r[3],'n,m)*av[m-1])/c);
 for(m=2,NBIG,
   my(c=subst(r[1],'n,m)); if(c==0, print("  b: leading coeff vanishes at n=",m); break);
   bv[m+1] = -(subst(r[2],'n,m)*bv[m] + subst(r[3],'n,m)*bv[m-1])/c);
 \\ integrality of d_n^k b_n
 my(dn=1, k2=1, k1=1, kfail=0);
 for(n=1,min(NBIG,400), dn=lcm(dn,n);
   if(denominator(dn^2*bv[n+1])!=1, k2=0);
   if(denominator(dn*bv[n+1])!=1, k1=0));
 print("  d_n^2 b_n in Z for n<=",min(NBIG,400),": ",k2,"    d_n b_n in Z: ",k1," (k=2 sharp iff 1/0)");
 \\ Casoratian
 my(W=vector(20,i,av[i+1]*bv[i+2]-av[i+2]*bv[i+1]));
 print("  Casoratian W_0..W_6 = ",vector(7,j,W[j]));
 \\ limit
 my(xi=bv[NBIG+1]/av[NBIG+1]*1.0, xi2=bv[NBIG]/av[NBIG]*1.0);
 print("  xi ~ ",xi,"   (prev ",xi2,",  |diff| ",abs(xi-xi2),")");
 my(l1=0,l2=0); my(rr2=rts);
 my(mods=vector(#rr2,i,abs(rr2[i]))); 
 print("  |roots| = ",mods);
 my(mx=vecmax(mods), mn=vecmin(mods));
 print("  lambda_1=",mx,"  lambda_2=",mn,"   score = log(1/lambda_2)-2 = ",log(1/mn)-2,
       "   budget = log(lambda_1)-2 = ",log(mx)-2);
 [lam, mr, rts, av, bv]};
