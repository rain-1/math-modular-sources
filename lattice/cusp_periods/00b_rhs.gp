/* 00b_rhs.gp -- inhomogeneity R_n = (n+1)^2 B_{n+1} - (11n^2+11n+3)B_n - n^2 B_{n-1}
   for the four Gamma_1(5) weight-3 Eisenstein directions.                    */
default(realprecision, 60);
NT = 120;
frob(a,b,c,nn) = {my(A=vector(nn+1)); A[1]=1; A[2]=b; for(n=1,nn-1, A[n+2] = ((a*n^2+a*n+b)*A[n+1] - c*n^2*A[n])/(n+1)^2); A;}
build2(a,b,c,nn) = {my(A=frob(a,b,c,nn), y0=Ser(A,'t,nn+1), P, K, g, qs, tq, Kq, Pq, Fq); P = 1 - a*'t + c*'t^2 + O('t^(nn+1)); K = P*y0^2; g = intformal((1/K - 1)/'t); qs = 't*exp(g); tq = serreverse(qs); Kq = ('t*deriv(tq,'t))/tq; Pq = 1 - a*tq + c*tq^2; Fq = sqrt(Kq/Pq); [tq, Fq];}
BF = build2(11,3,-1,NT);
xq = BF[1]; FF = BF[2];
peelx(S, xs, na) = {my(av=vector(na+1), G=S, xp=1+O('t^(NT+1))); for(n=0,na, my(cc=polcoeff(G,n)); av[n+1]=cc; if(cc!=0, G = G - cc*xp); xp = xp*xs); av;}
re4 = [1,0,0,-1,0]; im4 = [0,1,-1,0,0];
vv(v,d) = v[((d-1)%5)+1];
cout(v,m) = sumdiv(m, d, vv(v,d)*d^2);
cinn(v,m) = sumdiv(m, d, vv(v,m/d)*d^2);
eich(cv) = sum(m=1,NT, cv[m]/m^2*'t^m) + O('t^(NT+1));
NN = NT-2;
rhs(nm, cv) = {my(Th=eich(cv), B=peelx(FF*Th, xq, NN), R=vector(12)); print("== ", nm); for(n=0,11, my(bm=if(n==0,0,B[n]), b0=B[n+1], b1=B[n+2]); R[n+1] = (n+1)^2*b1 - (11*n^2+11*n+3)*b0 - n^2*bm); print("   R_0..R_11 = ", R); print("   ratios R_{n+1}/R_n = ", vector(9,i, if(R[i+1]!=0, 1.0*R[i+2]/R[i+1], "-")));}
rhs("R1 (outer 2Re psi4)", vector(NT,m,cout([2,0,0,-2,0],m)));
rhs("R2 (outer -2Im psi4)", vector(NT,m,cout([0,-2,2,0,0],m)));
rhs("Phi_D", vector(NT,m,cout([1,-2,2,-1,0],m)));
ph5 = (11+5*sqrt(5))/2;
rhs("R3 (inner 2Re)", vector(NT,m,cinn([2,0,0,-2,0],m)));
rhs("R4 (inner -2Im)", vector(NT,m,cinn([0,-2,2,0,0],m)));
print("t1 = ",(11-5*sqrt(5))/2*1.0,"  t2 = ",-(11+5*sqrt(5))/2*1.0, "  1/t2 = ", -2/(11+5*sqrt(5)));
quit;
