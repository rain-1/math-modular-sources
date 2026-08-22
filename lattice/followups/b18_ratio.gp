/* b18_ratio.gp -- allow xi to be a RATIO: test (a+b sqrt17 + c pi^2 + ...) * xi in the
   Gamma-class span, at both places.  Also a clean good-prime tower measurement.   */
default(parisize, 6000000000);
default(realprecision, 1100);
xi = eval(Str(read("b14_xi_arch_1200.txt")))*1.0;
s17=sqrt(17.); zm=(349-85*s17)/131072; Lam=log(abs(1/zm));
P2=Pi^2; Z3=zeta(3);
tv(nm,v,hb)={my(r=lindep(v,1000));
  if(type(r)=="t_COL" && r[1]!=0 && vecmax(abs(r))<hb, print("  *** ",nm," : ",r~), print("  ",nm," : none"));}
G = [1, Lam, Lam^2, Lam^3, P2, P2*Lam, Z3];
print("=== archimedean, xi with an algebraic/period denominator ===");
tv("[xi, s17*xi] + G",           concat([xi, s17*xi], G), 10^30);
tv("[xi, s17*xi] + G + s17*G",   concat([xi, s17*xi], concat(G, s17*G)), 10^20);
tv("[xi, P2*xi] + G",            concat([xi, P2*xi], G), 10^30);
tv("[xi, s17*xi, P2*xi, P2*s17*xi] + G", concat([xi,s17*xi,P2*xi,P2*s17*xi], G), 10^25);
tv("[xi,s17*xi] + {1,s17,Z3,s17*Z3,P2,s17*P2}", concat([xi,s17*xi],[1,s17,Z3,s17*Z3,P2,s17*P2]), 10^40);
tv("[xi] + {1,s17,Z3,s17*Z3}",   concat([xi],[1,s17,Z3,s17*Z3]), 10^60);
tv("[xi] + {1,Lam,Z3,P2*Lam,Lam^3}", concat([xi],[1,Lam,Z3,P2*Lam,Lam^3]), 10^50);
print("\n=== good-prime tower measurement (clean) ===");
read("../mum_survey/apery.gp"); read("../mum_survey/ops.gp");
OP=0; for(i=1,#OPS, if(OPS[i][1]=="207", OP=OPS[i]));
NNt=6000; prt=aperyPair(OP[4],NNt); At=prt[1]; Bt=prt[2];
{ for(j=1,3, my(p=[3,5,7][j]);
    for(a=1,3,
      my(smax=0); while(a*p^(smax+1)<=NNt, smax++);
      print("  p=",p,"  a=",a,"  smax=",smax);
      for(s=0,smax,
        my(n=a*p^s, t=Bt[n+1]/At[n+1]);
        printf("    s=%d n=%-5d v_p(A)=%-3d v_p(B)=%-4d v_p(t)=%-5d",
               s, n, valuation(At[n+1],p), valuation(Bt[n+1],p), valuation(t,p));
        if(s>0, my(rr = t/(Bt[a*p^(s-1)+1]/At[a*p^(s-1)+1]));
           printf("  v_p(ratio)=%-4d  ratio*p^3 = %s", valuation(rr,p), Str(rr*p^3+O(p^5))));
        print()))); }
quit
