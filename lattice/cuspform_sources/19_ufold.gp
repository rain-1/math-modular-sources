default(realprecision, 60);
read("lib.gp");
read("hosts.gp");
NQ = 700;
for(i=1,#HOSTS, my(h=HOSTS[i], N=h[1], C=h[2], us, qc, uc); us=useries(h[4],h[5],NQ); qc=exp(-2*Pi/sqrt(N*1.0)); uc=subst(truncate(us),q,qc); print(h[6], "   u(i/sqrtN) = ", uc, "   1/sqrt(C) = ", 1.0/sqrt(C*1.0), "   u*sqrt(C) = ", uc*sqrt(C*1.0)));
quit;
