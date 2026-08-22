\\ The zeta(7) level-24 system in the *s* coordinate (s=(z/(1-z))^2), where it is a
\\ genuine scalar Sym^6:  Atilde(s) = sqrt(1-34s+s^2)*Apery(s)^3.
\\ Its 6th root h = (1-34s+s^2)^{1/12} * sqrt(Apery(s)) is a SECOND-ORDER row.
read("lattice/root_rows/lib.gp");
default(realprecision,260);
N=320;
ap=vector(N+1,i,my(m=i-1);sum(k=0,m,binomial(m,k)^2*binomial(m+k,k)^2));
Aser=sum(m=0,N,ap[m+1]*'s^m)+O('s^(N+1));
PP=1-34*'s+'s^2+O('s^(N+1));
At=Aser^3*sqrt(PP);
Av=vector(N+1,i,polcoeff(At,i-1));
print("parent Atilde_n (n=0..8) = ",vector(9,j,Av[j]));
print("parent integral: ",sum(i=1,N+1,denominator(Av[i])!=1)==0);
rr=rootrow(Av,6,N); lam=rr[1]; a=rr[2];
print("lambda = ",lam);
print("a_0..a_6 = ",vector(7,j,a[j]));
print("a_n integral to n=",N,": ",sum(i=1,N+1,denominator(a[i])!=1)==0);
mr=minrec(a,4,4);
print("minimal recurrence: order=",mr[1]," deg=",mr[2]);
print("   ",mr[3]);
\\ operator form: L = sum_j s^j q_j(theta),  q_j(m) = r_j(m+j)
r=mr[3];
qs=vector(5,j,subst(r[j],'n,'n+(j-1)));
print("operator L' = sum_j s^j q_j(theta):");
{for(j=1,5, print("   q_",j-1,"(th) = ",qs[j]));}
lead = sum(j=0,4, polcoeff(qs[j+1],2,'n)*'S^j);
print("leading coefficient A(s) = ",lead,"   factored: ",factor(lead));
\\ companion: L'(B) = s, b_0=0, b_1=1
NB=320;
b=vector(NB+1); b[1]=0; b[2]=1;
{for(m=2,NB, my(t=0); for(j=1,4, if(m-j>=0, t+=subst(qs[j+1],'n,m-j)*b[m-j+1]));
  b[m+1] = -t/m^2);}
print("b_1..b_6 = ",vector(6,j,b[j+1]));
\\ denominators
{my(dn=1,k1=1,k2=1,k3=1,fail2=0);
 for(n=1,300, dn=lcm(dn,n);
   if(denominator(dn*b[n+1])!=1,k1=0);
   if(denominator(dn^2*b[n+1])!=1, if(k2==1,print("   d_n^2 b_n fails first at n=",n," denom=",denominator(dn^2*b[n+1]))); k2=0);
   if(denominator(dn^3*b[n+1])!=1,k3=0));
 print("d_n b_n in Z: ",k1,"   d_n^2 b_n in Z: ",k2,"   d_n^3 b_n in Z: ",k3);}
\\ the limit
xi = b[N+1]/a[N+1]*1.0;
print("xi ~ ",xi);
print("xi (n-1) ~ ",b[N]/a[N]*1.0);
print("|difference| = ",abs(b[N+1]/a[N+1]-b[N]/a[N])*1.0);
\\ rates
{my(rts=charroots(mr[3]), mods=vector(#rts,i,abs(rts[i])));
 print("|char roots| = ",mods);
 my(l1=vecmax(mods), l2=vecmin(mods));
 print("lambda_1 = ",l1,"  lambda_2 = ",l2);
 print("score(k=2) = ",log(1/l2)-2,"   score(k=3) = ",log(1/l2)-3);
 print("budget(k=2) = ",log(l1)-2);}
\\ empirical decay of the linear form
{for(n=60,300,if(n%60==0, print("  n=",n,"  log|a_n xi - b_n|/n = ", log(abs(a[n+1]*xi-b[n+1]))/n)));}
write("lattice/root_rows/zeta7_s_xi.txt", xi);
\q
