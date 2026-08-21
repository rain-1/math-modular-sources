\\ Master census: six Zagier order-2 sporadics + six AZ order-3 sporadics.
\\ Exact integer arithmetic. gp -q master12.gp
default(realprecision,120);

L2chi3=lfun(-3,2);
L3chi3=lfun(-3,3);
basis=[1,zeta(2),zeta(3),zeta(5),Catalan,L2chi3,L3chi3,Pi^3,Pi^4,log(2)];
bnames=["1","zeta(2)","zeta(3)","zeta(5)","G","L(2,chi-3)","L(3,chi-3)","pi^3","pi^4","log2"];

\\ order-2 rows: (n+1)^2 u_{n+1} = (a n^2+a n+b) u_n - c n^2 u_{n-1}, a0=1,a1=b, b0=0,b1=1
\\ integer model: A_n = a_n (integer); Cc_n with b_n = Cc_n / (n!)^2  (Cc_0=0,Cc_1=1)
row2(a,b,c,N)={
  my(A=vector(N+1),Cc=vector(N+1));
  A[1]=1;A[2]=b;Cc[1]=0;Cc[2]=1;
  for(n=1,N-1,
    A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;
    Cc[n+2]=(a*n^2+a*n+b)*Cc[n+1]-c*n^4*Cc[n];
  );
  [A,Cc]};

\\ order-3 rows: (n+1)^3 u_{n+1} = (2n+1)(a n^2+a n+b) u_n - c n^3 u_{n-1}
row3(a,b,c,N)={
  my(A=vector(N+1),Cc=vector(N+1));
  A[1]=1;A[2]=b;Cc[1]=0;Cc[2]=1;
  for(n=1,N-1,
    A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
    Cc[n+2]=(2*n+1)*(a*n^2+a*n+b)*Cc[n+1]-c*n^6*Cc[n];
  );
  [A,Cc]};

bval(Cc,ord,n)=Cc[n+1]/(n!)^ord;

denomexp(A,Cc,ord,N)={
  my(kmax=0,d);
  for(n=2,N,
    d=lcm(vector(n,i,i));
    my(bn=bval(Cc,ord,n));
    my(k=0);
    while(denominator(bn*d^k)!=1 && k<20, k++);
    if(k>kmax,kmax=k);
  );
  kmax};

slopeatp(A,Cc,ord,N,p)={
  my(rprev=bval(Cc,ord,1)/A[2],vs=vector(0));
  for(n=2,N,
    my(rn=bval(Cc,ord,n)/A[n+1]);
    my(dv=rn-rprev);
    if(dv!=0, vs=concat(vs,valuation(dv,p)));
    rprev=rn;
  );
  vs};

maxvp(A,N,p)={my(m=0);for(n=1,N,if(A[n+1]!=0,my(v=valuation(A[n+1],p));if(v>m,m=v)));m};

identify1(x)={
  my(xr=precision(x*1.0,80));
  my(v=vector(#basis+1)); for(i=1,#basis,v[i]=basis[i]*1.0); v[#basis+1]=-xr;
  lindep(v)};

\\ identify using two different truncations N1<N2; require x-coeff nonzero and stable
identify(Cc,A,ord,N1,N2)={
  my(x1=bval(Cc,ord,N1)/A[N1+1]);
  my(x2=bval(Cc,ord,N2)/A[N2+1]);
  my(r1=identify1(x1), r2=identify1(x2));
  if(r1[#basis+1]==0 || r2[#basis+1]==0, return(["unidentified (trivial/zero relation)",r1,r2]));
  if(r1 != r2 && r1 != -r2, return(["unstable/unidentified",r1,r2]));
  r1};

report(label,a,b,c,ord,N,plist)={
  my(row = if(ord==2,row2(a,b,c,N),row3(a,b,c,N)));
  my(A=row[1],Cc=row[2]);
  print("== ",label," (a,b,c)=(",a,",",b,",",c,") ord=",ord," ==");
  \\ characteristic roots from x^2-(leading a coeff... for order2: x^2-a x + c=0? use n->infty coeff ratio
  my(pol = if(ord==2, x^2-a*x+c, x^2-2*a*x+c));
  my(rt=polroots(pol));
  print("  char roots (from x^2-",if(ord==2,a,2*a),"x+",c,"=0): ",rt[1],", ",rt[2]);
  my(k=denomexp(A,Cc,ord,N));
  print("  denom exponent k (n<=",N,") = ",k);
  for(pi=1,#plist,
    my(p=plist[pi]);
    my(vs=slopeatp(A,Cc,ord,N,p));
    my(last=if(#vs>0,vs[#vs],0));
    print("  slope at p=",p,": last few valuations = ",if(#vs>5,vs[#vs-4..#vs],vs)," predicted v_p(c)=",valuation(c,p));
  );
  my(mv2=maxvp(A,N,2),mv3=maxvp(A,N,3),mv5=maxvp(A,N,5));
  print("  max v_2(a_n)=",mv2," max v_3(a_n)=",mv3," max v_5(a_n)=",mv5);
  my(rt2=polroots(if(ord==2,x^2-a*x+c,x^2-2*a*x+c)));
  my(iscomplex = abs(imag(rt2[1]))>1e-20);
  if(iscomplex,
    print("  [complex conjugate roots -> no real archimedean limit; skipping lindep]"),
    my(lin=identify(Cc,A,ord,N-20,N));
    print("  lindep vs basis ",bnames," -> ",lin);
  );
  [A,Cc]
};

plist=[2,3,5,7,11,13];

print("###### ORDER-2 ZAGIER ROWS ######");
r_A=report("Zagier A",7,2,-8,2,300,plist);
r_B=report("Zagier B",9,3,27,2,300,plist);
r_C=report("Zagier C",10,3,9,2,300,plist);
r_D=report("Zagier D",11,3,-1,2,300,plist);
r_E=report("Zagier E",12,4,32,2,300,plist);
r_F=report("Zagier F",17,6,72,2,300,plist);

print("###### ORDER-3 AZ ROWS ######");
r_del=report("AZ delta(7,3,81)",7,3,81,3,220,plist);
r_ACnn=report("AZ(9,3,-27)",9,3,-27,3,220,plist);
r_Domb=report("Domb(10,4,64)",10,4,64,3,220,plist);
r_eta=report("AZ eta(11,5,125)",11,5,125,3,220,plist);
r_zz=report("AZ(12,4,16) [T]",12,4,16,3,220,plist);
r_Apery=report("Apery(17,5,1)",17,5,1,3,220,plist);
