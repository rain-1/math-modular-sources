/* Mirror maps t_X(q) for Zagier rows, from the recurrence alone.
   Row: (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1}, P = aa*(n^2+n)+bb, Q = dd*n^2.
   L = th^2 - t P(th) + t^2 Q(th+1).   Variable 'x' plays the role of t, then q. */
NT = 30;

{ rowu(aa,bb,dd,N) =
  my(u=vector(N+2));
  u[1]=1;
  for(n=0,N,
      my(Pn=aa*(n^2+n)+bb, Qn=dd*n^2, prv=if(n==0,0,u[n]));
      u[n+2] = (Pn*u[n+1] - Qn*prv)/(n+1)^2);
  vector(N+1,k,u[k]);
}

{ rowv(aa,bb,dd,N) =
  my(u=rowu(aa,bb,dd,N), v=vector(N+2));
  for(n=1,N,
      my(P1=aa*((n-1)^2+(n-1))+bb, Qm=dd*(n-1)^2);
      my(Pp=2*aa*(n-1)+aa, Qp=2*dd*(n-1));
      my(U=u[n+1], Um=u[n], Um2=if(n>=2,u[n-1],0));
      my(Vm=v[n], Vm2=if(n>=2,v[n-1],0));
      my(rhs = -(2*n*U - Pp*Um + Qp*Um2));
      v[n+1] = (rhs + P1*Vm - Qm*Vm2)/n^2);
  vector(N+1,k,v[k]);
}

{ mirror(aa,bb,dd,N) =
  my(u=rowu(aa,bb,dd,N), v=rowv(aa,bb,dd,N));
  my(y0=Ser(u,'x,N+1), h=Ser(v,'x,N+1));
  my(qq='x*exp(h/y0));
  serreverse(qq);
}

Arow=[7,2,-8]; Crow=[10,3,9]; Frow=[17,6,72];
tA = mirror(Arow[1],Arow[2],Arow[3],NT);
tC = mirror(Crow[1],Crow[2],Crow[3],NT);
tF = mirror(Frow[1],Frow[2],Frow[3],NT);
print("t_A = ", tA + O('x^9));
print("t_C = ", tC + O('x^9));
print("t_F = ", tF + O('x^9));
print();
print("A - C/(1-C)      : ", (tA - tC/(1-tC)) + O('x^41));
print("F + C/(1-9C)     : ", (tF + tC/(1-9*tC)) + O('x^41));
print("F - C/(1-9C)     : ", (tF - tC/(1-9*tC)) + O('x^10));
print("F - A/(1+8A)     : ", (tF - tA/(1+8*tA)) + O('x^10));
print("F + A/(1-8A)     : ", (tF + tA/(1-8*tA)) + O('x^10));
print();
print("/* F - w is even and F + w is odd, w := t_C/(1-9 t_C): hence w(-x) = -t_F(x). */");
w = tC/(1-9*tC);
print("t_F(x) + w(-x)   : ", tF + subst(w,'x,-'x));
print("  ==> the F-row is the (-1)^n twist of the third X_0(6) coordinate w.");
quit;
