/* verifies every closed form asserted in consolidation/ZUDILIN_2ADIC.md */
polpq(N)={my(P=vector(N+2),Q=vector(N+2));Q[1]=1;Q[2]=X^2-X+1;P[1]=0;P[2]=1;for(m=1,N, Q[m+2]=((2*m*(m+1)+1-X+X^2)*Q[m+1]-m^2*Q[m])/(m+1)^2; P[m+2]=((2*m*(m+1)+1-X+X^2)*P[m+1]-m^2*P[m])/(m+1)^2);[P,Q];}
read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/zud_row.gp");
NM=24; b=polpq(NM+2); PP=b[1]; QQ=b[2];
sh(f,c)=subst(f,X,X+c);
H(m)=sum(k=1,m,1/k);
PZ(n)=20*n^2-8*n+1;
QZ(n)=3520*n^6+5632*n^5+2064*n^4-384*n^3-156*n^2+16*n+7;
Wp(n)=X^2*(X-1)^2*(sh(PP[n+2],-1)*sh(QQ[n],1) - sh(PP[n],1)*sh(QQ[n+2],-1)) - (4*X-2)*sh(QQ[n+2],-1)*sh(QQ[n],1);
Wcl(n)=(X^6-(4*n+5)*X^5+10*(n+1)^2*X^4-(2*n+1)*(2*n+3)*(4*n+3)*X^3+(2*n+1)*(2*n+3)*(4*n^2+6*n+1)*X^2-2*n*(n+1)*(4*n^3+14*n^2+14*n+3)*X+2*n^2*(n+1)^2*(2*n+3))/(n^2*(n+1)^2);
Om(n)=X^2*(PP[n+1]*sh(QQ[n],1)+sh(PP[n],1)*QQ[n+1])-2*QQ[n+1]*sh(QQ[n],1);
qh(n)=subst(QQ[n+1],X,1/2-n);
chk(nm,f)={my(ok=1); for(n=1,NM, if(!f(n), ok=0)); print("  ",nm,"  ",if(ok,"OK","*** FAIL ***"));}
r=zudrow(NM+2); u=r[1]; v=r[2];
chk("W_n closed form as typed", n->Wp(n)==Wcl(n));
chk("W_n at x_n  = q(n)/[64 n^2 (n+1)^2]", n->subst(Wp(n),X,1/2-n)==QZ(n)/(64*n^2*(n+1)^2));
chk("Omega'_n at x_n = q(n)/[4 n^2 (n+1)^2 (2n-1)^2 (2n+1)^2]", n->subst(Wp(n),X,1/2-n)/((1/2-n)^2*(-1/2-n)^2)==QZ(n)/(4*n^2*(n+1)^2*(2*n-1)^2*(2*n+1)^2));
chk("W_n at 0  = 2(2n+3)", n->subst(Wp(n),X,0)==2*(2*n+3));
chk("W_n at 1  = -2(2n-1)", n->subst(Wp(n),X,1)==-2*(2*n-1));
chk("W_n' at 0", n->subst(deriv(Wp(n),X),X,0)==-(2*n+3)*(4+2/n+2/(n+1)));
chk("W_n' at 1", n->subst(deriv(Wp(n),X),X,1)==(2*n-1)*(-4+2/n+2/(n+1)));
chk("W_n at 1/2 = [8n^2+8n+7]/[64n^2(n+1)^2]", n->subst(Wp(n),X,1/2)==(8*n^2+8*n+7)/(64*n^2*(n+1)^2));
chk("Lemma A: x^2 Omega_n = -[x^2/n^2 - 2x/n + 2]", n->Om(n)==-(X^2/n^2-2*X/n+2));
chk("Lemma A at x_n = -p(n)/[n^2(2n-1)^2]", n->subst(Om(n),X,1/2-n)/(1/2-n)^2 == -PZ(n)/(n^2*(2*n-1)^2));
chk("q_m at 0 = 1", m->subst(QQ[m+1],X,0)==1);
chk("q_m at -1 = 2m+1", m->subst(QQ[m+1],X,-1)==2*m+1);
chk("q_m at 2 = 2m+1", m->subst(QQ[m+1],X,2)==2*m+1);
chk("q_m' at 0 = -H_m", m->subst(deriv(QQ[m+1],X),X,0)==-H(m));
chk("q_m' at -1 = -[2m+1]H_m", m->subst(deriv(QQ[m+1],X),X,-1)==-(2*m+1)*H(m));
chk("Casoratian p_n q_n-1 - q_n p_n-1 = 1/n^2", n->PP[n+1]*QQ[n]-QQ[n+1]*PP[n]==1/n^2);
chk("derived 3-term relation for qhat", n-> PZ(n)/(n^2*(2*n-1)^2)*qh(n+1) == QZ(n)/(4*n^2*(n+1)^2*(2*n-1)^2*(2*n+1)^2)*qh(n) + PZ(n+1)/((n+1)^2*(2*n+1)^2)*qh(n-1));
chk("qhat satisfies Zudilin recursion", n-> (2*n+1)^2*(2*n+2)^2*PZ(n)*qh(n+1) == QZ(n)*qh(n)+(2*n-1)^2*(2*n)^2*PZ(n+1)*qh(n-1));
chk("Zudilin Casoratian = (-1)^(m-1) p(m)/[8 m^2 (2m-1)^2]", m-> v[m+1]*u[m]-v[m]*u[m+1] == (-1)^(m-1)*PZ(m)/(8*m^2*(2*m-1)^2));
quit;
