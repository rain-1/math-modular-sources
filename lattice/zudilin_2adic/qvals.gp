polpq(N)={my(P=vector(N+2),Q=vector(N+2));Q[1]=1;Q[2]=X^2-X+1;P[1]=0;P[2]=1;
 for(m=1,N, Q[m+2]=((2*m*(m+1)+1-X+X^2)*Q[m+1]-m^2*Q[m])/(m+1)^2;
            P[m+2]=((2*m*(m+1)+1-X+X^2)*P[m+1]-m^2*P[m])/(m+1)^2);[P,Q];}
N=14;b=polpq(N);P=b[1];Q=b[2];
H(m)=sum(k=1,m,1/k);
print("m  q_m(0) q_m(-1) q_m(-2) q_m'(0)+H(m)  q_m'(-1)");
{for(m=0,N, my(q=Q[m+1], d=deriv(q,X));
  print("  ",m,"  ",subst(q,X,0)," ",subst(q,X,-1)," ",subst(q,X,-2),"  ",subst(d,X,0)+H(m),"   ",subst(d,X,-1)));}
print("");
print("guess q_m'(-1) = -(2m+1)H(m) + something:");
{for(m=0,N, my(d=deriv(Q[m+1],X), v=subst(d,X,-1)); print("  m=",m,"  q'(-1)=",v,"  v+(2m+1)*H(m)=",v+(2*m+1)*H(m), "   v+(2m+1)*H(m)-2*m=",v+(2*m+1)*H(m)-2*m));}
quit;
