\\ census_scores.gp -- p-adic irrationality score S_p for every census cell.
\\ Conventions (see consolidation/PADIC_IRRATIONALITY_CENSUS.md sec.1):
\\   A_n, B_n : the row's two solutions, B_n/A_n -> xi_p in Q_p.
\\   kap  = -lim v_p(A_n)/n   (rate of p-power DENOMINATORS; 0 for integral rows)
\\   sg   =  lim v_p(xi_p - B_n/A_n)/n   (quality slope; = v_p(c) + 2*kap)
\\   kk   =  lim log(clearing factor)/n  (prime-to-p lcm cost, nats per index)
\\   lam1 =  lim log max(|A_n|,|B_n|)/n
\\   S_p  =  sg*log p - kk - kap*log p - lam1      (irrational if S_p > 0)
\\   th_p =  sg*log p / (kk + kap*log p + lam1)    (Calegari's theta; th>1 <=> S>0)
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/census_scores.gp
default(parisizemax,"8G");
default(realprecision,60);

NORD2 = 500;
NORD3 = 420;

\\ order-2 Zagier: (n+1)^2 u_{n+1} = (a n^2 + a n + b) u_n - c n^2 u_{n-1}
{ row2(a,b,c,N) =
  my(A=vector(N+1),Cc=vector(N+1));
  A[1]=1; A[2]=b; Cc[1]=0; Cc[2]=1;
  for(n=1,N-1,
    A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;
    Cc[n+2]=(a*n^2+a*n+b)*Cc[n+1]-c*n^4*Cc[n];
  );
  [A,Cc,2];
}

\\ order-3 AZ: (n+1)^3 u_{n+1} = (2n+1)(a n^2+a n+b) u_n - c n^3 u_{n-1}
{ row3(a,b,c,N) =
  my(A=vector(N+1),Cc=vector(N+1));
  A[1]=1; A[2]=b; Cc[1]=0; Cc[2]=1;
  for(n=1,N-1,
    A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
    Cc[n+2]=(2*n+1)*(a*n^2+a*n+b)*Cc[n+1]-c*n^6*Cc[n];
  );
  [A,Cc,3];
}

\\ (Cooper rows: see cooper_scores.gp -- their recurrence is not in Zagier/AZ form)

bval(R,n) = R[2][n+1]/(n!)^R[3];

\\ smallest integer kk with lcm(1..n)^kk * B_n in Z, for all n <= N
{ kexp(R,N) =
  my(kmax=0,d=1);
  for(n=1,N, d=lcm(d,n);
    my(x=bval(R,n), kk=0);
    while(denominator(x*d^kk)!=1 && kk<15, kk++);
    if(kk>kmax, kmax=kk));
  kmax;
}

\\ measured kappa_p and its stability
{ kapmeas(R,pp,N) =
  my(v=vector(4)); for(j=1,4, my(n=N*j\4); v[j]=-valuation(R[1][n+1],pp)*1./n); v;
}

\\ measured quality slope sigma_p
{ sigslope(R,pp,N,Nmeas) =
  my(xi=bval(R,N)/R[1][N+1], n1=Nmeas\2, n2=Nmeas);
  my(v1=valuation(xi-bval(R,n1)/R[1][n1+1],pp), v2=valuation(xi-bval(R,n2)/R[1][n2+1],pp));
  (v2-v1)*1./(n2-n1);
}

{ archmeas(R,N) =
  my(out=vector(3));
  for(j=1,3, my(n=N-(3-j)*(N\6)); out[j]=log(max(abs(R[1][n+1]*1.),abs(bval(R,n)*1.)))/n);
  out;
}

\\ ------------------------------------------------------------------ the census
\\ name, order, (a,b,c), primes to test, xi_p identification
DATA2 = [ ["A",7,2,-8], ["B",9,3,27], ["C",10,3,9], ["D",11,3,-1], ["E",12,4,32], ["F",17,6,72] ];
DATA3 = [ ["delta",7,3,81], ["zeta",9,3,-27], ["Domb",10,4,64], ["eta",11,5,125], ["T=eps",12,4,16], ["Apery",17,5,1] ];
PR = [2,3,5,7];

print("row  p  kappa_p  sigma_p(meas)  v_p(c)  k  log_lambda1  S_p  theta_p");
{
for(i=1,#DATA2,
  my(d=DATA2[i], nm=d[1], a=d[2], b=d[3], c=d[4], N=NORD2);
  my(R=row2(a,b,c,N));
  my(kk=kexp(R,200));
  my(disc=a^2-4*c, lam1=if(disc>=0, (abs(a)+sqrt(disc*1.))/2, sqrt(abs(c*1.))));
  my(am=archmeas(R,N));
  for(j=1,#PR, my(pp=PR[j]);
    if(valuation(c,pp)>0,
      my(kap=kapmeas(R,pp,N), sg=sigslope(R,pp,N,400));
      my(kapv=kap[4]);
      my(Sp = sg*log(pp) - kk - kapv*log(pp) - log(lam1));
      my(th = sg*log(pp) / (kk + kapv*log(pp) + log(lam1)));
      print(nm,"  p=",pp,"  kappa~",kap[4],"  sigma_meas=",sg,"  v_p(c)=",valuation(c,pp),
            "  k=",kk,"  log_l1=",log(lam1),"  arch_meas=",am[3],"  S=",Sp,"  theta=",th);
    ));
);
for(i=1,#DATA3,
  my(d=DATA3[i], nm=d[1], a=d[2], b=d[3], c=d[4], N=NORD3);
  my(R=row3(a,b,c,N));
  my(kk=kexp(R,200));
  my(disc=4*a^2-4*c, lam1=if(disc>=0, (abs(2*a)+sqrt(disc*1.))/2, sqrt(abs(c*1.))));
  my(am=archmeas(R,N));
  for(j=1,#PR, my(pp=PR[j]);
    if(valuation(c,pp)>0,
      my(kap=kapmeas(R,pp,N), sg=sigslope(R,pp,N,400));
      my(kapv=kap[4]);
      my(Sp = sg*log(pp) - kk - kapv*log(pp) - log(lam1));
      my(th = sg*log(pp) / (kk + kapv*log(pp) + log(lam1)));
      print(nm,"  p=",pp,"  kappa~",kap[4],"  sigma_meas=",sg,"  v_p(c)=",valuation(c,pp),
            "  k=",kk,"  log_l1=",log(lam1),"  arch_meas=",am[3],"  S=",Sp,"  theta=",th);
    ));
);
}

print("\n(Cooper s7,s10,s18: run cooper_scores.gp)");
quit;
