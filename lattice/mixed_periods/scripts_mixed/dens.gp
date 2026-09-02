N=160;
lcmv(n)=my(l=1);for(k=1,n,l=lcm(l,k));l;
intg(s)=my(v=Vec(s+O(x^(N+1))));sum(k=1,#v, v[k]*x^k/k)+O(x^(N+1)); \\ integral from 0
den(s)=my(v=Vec(s+O(x^(N+1))));vector(#v,k,denominator(v[k]));
\\ report: for each n, is den | L(n)^a L(n/2)^b ? report exponents needed minimal pattern
report(name,s)={my(d=den(s),ok1=1,ok11=1,ok2=1,bad=[]);
 for(n=1,N-1, my(L=lcmv(n),L2=lcmv(n\2), dn=d[n+1]);
   if(dn%1==0, );
   if(L%dn!=0, ok1=0);
   if((L*L2)%dn!=0, ok11=0; bad=concat(bad,[n]));
   if((L*L)%dn!=0, ok2=0));
 print(name,": [1..n]:",ok1,"  [1..n][1..n/2]:",ok11,"  [1..n]^2:",ok2, "  first fails of [n][n/2]: ",if(#bad,vector(min(6,#bad),i,bad[i]),"none"));
}
host(name,HA)={
 print("=== host ",name);
 my(lg=log(1-x+O(x^(N+2))));
 HB=HA*intg(HA/(1-x)); HC=HA*intg(HA*lg/x); HD=HA*intg(HA*lg/(1-x)); HE=HA*intg(HA*lg^2/x); HF=HA*intg(HA/x*lg)*0;
 HG=HA*intg(HA*lg^2/(1-x));
 HH=HA*intg(HA*intg(HA/(1-x))/HA/x);
 report("HA",HA);report("HB",HB);report("HC",HC);report("HD",HD);report("HE log^2/t",HE);report("HG log^2/(1-t)",HG);
 \\ also plain integrals without HA prefactor
 report("int HA/(1-t)",intg(HA/(1-x)));report("int HA log/t",intg(HA*lg/x));report("int HA log/(1-t)",intg(HA*lg/(1-x)));
}
host("1/sqrt(1-4x)",1/sqrt(1-4*x+O(x^(N+2))));
host("1/sqrt(1-8x)",1/sqrt(1-8*x+O(x^(N+2))));
host("1/sqrt(1-12x)",1/sqrt(1-12*x+O(x^(N+2))));
host("1/sqrt((1-x)(1-9x))",1/sqrt((1-x)*(1-9*x)+O(x^(N+2))));
host("1/sqrt((1-x)(1-25x))",1/sqrt((1-x)*(1-25*x)+O(x^(N+2))));
host("1/sqrt(1-6x+x^2) Delannoy",1/sqrt(1-6*x+x^2+O(x^(N+2))));
