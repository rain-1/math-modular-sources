default(realprecision, 120);
read("lib.gp");
read("hosts.gp");
NQ = 700;
for(i=1,#HOSTS, my(h=HOSTS[i], N=h[1], C=h[2], B=h[3], us, xs, qc, xc); us=useries(h[4],h[5],NQ); xs=us/(1+B*us+C*us^2); qc=exp(-2*Pi/sqrt(N*1.0)); xc=subst(truncate(xs),q,qc); print(h[6], "   x(i/sqrtN) = ", xc, "    1/lam1 = ", 1.0/(B+2*sqrt(C*1.0))));
quit;
