read("lib.gp");
iterrec(P, a0, a1, nb) = my(v=vector(nb+2)); v[1]=a0; v[2]=a1; for(n=1, nb, my(p0=subst(P[1],x,n), p1=subst(P[2],x,n), p2=subst(P[3],x,n)); if(p2==0, break); v[n+2] = -(p0*v[n]+p1*v[n+1])/p2); v;
