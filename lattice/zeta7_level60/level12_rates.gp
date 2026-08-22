\\ Level-12 zeta(7) parent: growth rates and minimal recurrence.
default(parisizemax, 8000000000);
default(parisize, 4000000000);
default(realprecision, 80);
N = 460;
et(k) = eta(q^k + O(q^(N+2)));
d12=[1,2,3,4,6,12]; c12=[1,-572,11583,-36608,46332,-20736];
{g = vector(N, m, sum(i=1,6, if(m%d12[i]==0, c12[i]*sigma(m\d12[i],7), 0)));}
Phi = sum(m=1,N, g[m]*q^m) + O(q^(N+1));
Psi = sum(m=1,N, (g[m]/m^7)*q^m) + O(q^(N+1));
xi  = (209/1728)*zeta(7);
x = q*et(4)^2*et(12)^2/(et(1)^2*et(3)^2);
t = x/(16*x^2+2*x+1);
Dt = q*deriv(t,q); Qt = serreverse(t);
A  = subst(Phi/Dt, q, Qt);
B  = A*subst(Psi, q, Qt);
M  = N-12;
An = vector(M+1,i,polcoeff(A,i-1));
Bn = vector(M+1,i,polcoeff(B,i-1));
print("rows built to n=", M);
print("A_n integral? ", vecmax(vector(M+1,i,denominator(An[i])))==1);
{for(i=1,9, my(n=[50,100,150,200,250,300,350,400,M][i]);
  print("n=",n,"  log|A_n|/n = ", log(abs(1.0*An[n+1]))/n,
        "   log|B_n - xi A_n|/n = ", log(abs(1.0*Bn[n+1]-xi*An[n+1]))/n,
        "   |B/A - xi|^(1/n) = ", abs(1.0*Bn[n+1]/An[n+1]-xi)^(1.0/n)));}
write("lattice/zeta7_level60/level12_rows.txt", Str("An = ", An, ";"));
write("lattice/zeta7_level60/level12_rows.txt", Str("Bn = ", Bn, ";"));
\\ ---- mod-p scan for the minimal recurrence of A_n ----
PP = 2305843009213693951;
Am = vector(M+1, i, Mod(An[i], PP));
mk(V,r,D,nlo,nhi) = matconcat(vector(nhi-nlo+1, i, my(n=nlo+i-1); vector((r+1)*(D+1), c, my(j=(c-1)\(D+1), k=(c-1)%(D+1)); n^k*V[n-j+1]))~);
print("mod-p kernel dimensions (order r, degree D):");
{for(rr=6,20, my(s="");
  for(DD=3,14, my(nv=(rr+1)*(DD+1), nhi=min(M, rr+1+nv+40));
    if(nhi-rr-1 >= nv+25, s=Str(s,"  D=",DD,":",#matker(mk(Am,rr,DD,rr+1,nhi))), s=Str(s,"  D=",DD,":-")));
  print("  r=",rr,s));}
print("DONE");
\q
