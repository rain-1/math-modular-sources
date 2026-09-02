/* 00_disc.gp -- settle the R1 / R2 discrepancy on Gamma_1(5) (Zagier D host). */
default(realprecision, 80);
NT = 300;
frob(a,b,c,nn) = {my(A=vector(nn+1)); A[1]=1; A[2]=b; for(n=1,nn-1, A[n+2] = ((a*n^2+a*n+b)*A[n+1] - c*n^2*A[n])/(n+1)^2); A;}
build2(a,b,c,nn) = {my(A=frob(a,b,c,nn), y0=Ser(A,'t,nn+1), P, K, g, qs, tq, Kq, Pq, Fq); P = 1 - a*'t + c*'t^2 + O('t^(nn+1)); K = P*y0^2; g = intformal((1/K - 1)/'t); qs = 't*exp(g); tq = serreverse(qs); Kq = ('t*deriv(tq,'t))/tq; Pq = 1 - a*tq + c*tq^2; Fq = sqrt(Kq/Pq); [tq, Fq];}
BF = build2(11,3,-1,NT);
xq = BF[1]; FF = BF[2];
print("x = ", xq + O('t^8));
print("F = ", FF + O('t^8));
peelx(S, xs, na) = {my(av=vector(na+1), G=S, xp=1+O('t^(NT+1))); for(n=0,na, my(cc=polcoeff(G,n)); av[n+1]=cc; if(cc!=0, G = G - cc*xp); xp = xp*xs); av;}
re4 = [1,0,0,-1,0]; im4 = [0,1,-1,0,0];
vv(v,d) = v[((d-1)%5)+1];
cout(v,m) = sumdiv(m, d, vv(v,d)*d^2);
cinn(v,m) = sumdiv(m, d, vv(v,m/d)*d^2);
eich(cv) = sum(m=1,NT, cv[m]/m^2*'t^m) + O('t^(NT+1));
phiD = [1,-2,2,-1,0];
R1v = [2,0,0,-2,0];
R2v = [0,-2,2,0,0];
ph5 = (11+5*sqrt(5))/2;
vnew = vector(5,j, 2*re4[j] - 2*ph5*im4[j]);
AA = peelx(FF, xq, NT-2);
print("A_n = ", vector(8,i,AA[i]));
NN = NT-2;
report(nm, cv, tgt) = {my(Th=eich(cv), B=peelx(FF*Th, xq, NN)); print("== ", nm); print("   B_1..5 = ", vector(5,i,B[i+1])); print("   B_n/A_n : n=",NN," ",B[NN+1]/AA[NN+1]*1.0); print("            n=",NN-50," ",B[NN-49]/AA[NN-49]*1.0); print("            n=",NN-150," ",B[NN-149]/AA[NN-149]*1.0); print("   target  = ", tgt); my(dn=B[NN+1]/AA[NN+1]*1.0-tgt, dm=B[NN-49]/AA[NN-49]*1.0-tgt); print("   defect at n=",NN,": ",dn,"   at n=",NN-50,": ",dm); print("   ratio of defects = ", dn/dm); B;}
BR1 = report("R1 (outer 2Re psi4)", vector(NT,m,cout(R1v,m)), 6*zeta(2)/5);
BR2 = report("R2 (outer -2Im psi4)", vector(NT,m,cout(R2v,m)), -2*zeta(2)/5);
BD  = report("Phi_D (outer w)", vector(NT,m,cout(phiD,m)), zeta(2)/5);
BN  = report("Phi_new (inner v)", vector(NT,m,cinn(vnew,m)), 0.6556341888406567663309814138723994024111);
quit;
