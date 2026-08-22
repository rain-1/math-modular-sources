\\ nf_herfurtner_uniformisation.gp
\\ Exact modular uniformisation of the two NEW Herfurtner rows of
\\ HERFURTNER_CLASSIFICATION.md 6.2.   (Fable, 2026-08-22)
\\
\\ Row #30  (n+1)^2 u_{n+1} = (117n^2+78n+21)u_n - 441(3n-1)^2 u_{n-1}, u_1=21
\\   t = -(1/9)*(eta(7 tau)/eta(tau))^4   (Gamma_0(7) hauptmodul, simple zero at ioo)
\\   A(t)*(1-117t+3969t^2)^(1/3) = theta_{-7} = 1 + 2 sum_n (sum_{d|n} chi_{-7}(d)) q^n
\\   equivalently  q dt/dq / t = theta_{-7}^2   (weight 2 on Gamma_0(7), chi_{-7}^2 = 1)
\\
\\ Row #45  (n+1)^2 u_{n+1} = (72n^2+36n+6)u_n - 108(4n-1)(4n-3)u_{n-1}, u_1=6
\\   j = (24t-1)^3/(8t^3),  i.e.  t = -1/(2(j^(1/3)-12)),  j^(1/3) = E_4/eta^8
\\   3*(q dt/dq)/t = A(t)^2 * (1-72t+1728t^2)^(1/2)     (w = q^(1/3), width-3 cusp)
default(parisizemax,6000000000);
default(seriesprecision,120);

print("==== Row #30, Gamma_0(7) ====");
{
my(NQ=90, av);
av=vector(NQ+3); av[1]=1; av[2]=21;
for(n=1,NQ+1, av[n+2]=((117*n^2+78*n+21)*av[n+1]-441*(3*n-1)^2*av[n])/(n+1)^2);
my(Aser=sum(n=0,NQ+1, av[n+1]*t^n)+O(t^(NQ+2)));
my(Rser=1-117*t+3969*t^2+O(t^(NQ+2)));
my(e1=eta(q+O(q^(NQ+3))), e7);
e7 = subst(e1,q,q^7)+O(q^(NQ+3));
my(tq = -q*(e7/e1)^4/9 + O(q^(NQ+3)));
my(Fq = subst(Aser*Rser^(1/3), t, tq));
my(E1 = 1+O(q^(NQ+1)));
for(n=1,NQ, my(s=0); fordiv(n,d, s+=kronecker(-7,d)); E1 += 2*s*q^n);
print("  a_n (n<=8)          : ",vector(9,j,av[j]));
print("  A*R^(1/3) - theta_-7 = ",(Fq-E1)+O(q^(NQ+1)));
print("  (q dt/dq)/t - theta^2= ",(q*deriv(tq,q)/tq - E1^2)+O(q^(NQ+1)));
}

print();
print("==== Row #45, level-3 index-3 group ====");
{
my(NT=45, av);
av=vector(NT+3); av[1]=1; av[2]=6;
for(n=1,NT+1, av[n+2]=((72*n^2+36*n+6)*av[n+1]-108*(4*n-1)*(4*n-3)*av[n])/(n+1)^2);
my(Aser=sum(n=0,NT+1, av[n+1]*t^n)+O(t^(NT+2)));
my(Rser=1-72*t+1728*t^2+O(t^(NT+2)));
my(e1=eta(w^3+O(w^(3*NT+6))));
my(E4=1+240*sum(n=1,NT+2, sigma(n,3)*w^(3*n))+O(w^(3*NT+6)));
my(jcr=E4/(w*e1^8));                 \\ j^(1/3)
my(tw=-1/(2*(jcr-12)));
my(Gw=subst(Aser^2*Rser^(1/2),t,tw));
my(Dt=(w/3)*deriv(tw,w));
print("  a_n (n<=8)          : ",vector(9,j,av[j]));
print("  t(w), w=q^(1/3)     : ",tw+O(w^8));
print("  3(q dt/dq)/t - A^2 sqrt(R) = ",3*Dt/tw - Gw + O(w^NT));
}
\q
