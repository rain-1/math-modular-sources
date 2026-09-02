/* 17_farperiods.gp -- local monodromy of A, B_D, B_new, B'_new at the two finite
   singular points t1 = phi^-5 and t2 = -phi^5, by numerical analytic continuation.
   For a solution y with no inhomogeneity pole at t, monodromy around t gives
   y -> y + 2 pi i * l(y) * u  with u the analytic local solution; the "period at t"
   is l(y)/l(A), i.e. the unique lambda with y - lambda A regular at t.          */
default(parisizemax, 20000000000);
default(realprecision, 200);
NC = 700; MT = 420;
s5 = sqrt(5); ph5 = (11+5*s5)/2; phm5 = (5*s5-11)/2; t1 = phm5; t2 = -ph5;
/* exact Taylor coefficients at x=0 */
ff = vector(NC+2); ff[1]=1; ff[2]=11; for(m=2,NC+0, ff[m+1]=11*ff[m]+ff[m-1]);
rD = vector(NC+1, i, if(i==1, 1, 0));
r3 = vector(NC+1, i, 2*ff[i]);
r4 = vector(NC+1, i, if(i==1, 0, -2*ff[i-1]));
mk(rr) = my(b=vector(NC+1)); b[1]=0; for(n=1,NC, b[n+1] = ((11*n^2-11*n+3)*b[n] + if(n>=2,(n-1)^2*b[n-1],0) + rr[n])/n^2); b;
Av = vector(NC+1); Av[1]=1; for(n=1,NC, Av[n+1] = ((11*n^2-11*n+3)*Av[n] + if(n>=2,(n-1)^2*Av[n-1],0))/n^2);
BDv = mk(rD); B3v = mk(r3); B4v = mk(r4);
print("check A_0..A_6 : ", vector(7,i,Av[i]));
print("check BD_1..BD_4: ", vector(4,i,BDv[i+1]));
print("check B3_1..B3_4: ", vector(4,i,B3v[i+1]));
print("check B4_1..B4_4: ", vector(4,i,B4v[i+1]));
/* numeric coefficient vectors */
Bnw  = vector(NC+1, i, B3v[i] + ph5*B4v[i]);
Bnwp = vector(NC+1, i, B3v[i] - phm5*B4v[i]);
ev(cf, x)  = my(s=0); forstep(k=NC, 0, -1, s = s*x + cf[k+1]); s;
evd(cf, x) = my(s=0); forstep(k=NC, 1, -1, s = s*x + k*cf[k+1]); s;
/* local RHS series generators: R(x) = c/(e0+e1*(x-x0)) or constant */
rser0(x0, M) = vector(M+1, i, 0);
rserC(x0, M, c) = vector(M+1, i, if(i==1, c, 0));
rserP(x0, M, c, a) = my(e0 = 1-a*x0); my(v=vector(M+1)); my(t=c/e0); for(i=0,M, v[i+1]=t; t = t*a/e0); v;
/* one Taylor step of the ODE  p y'' + q y' + r y = R,  p = x-11x^2-x^3 etc. */
stp(y0, y1, x0, h, Rv) = my(p0=x0-11*x0^2-x0^3, p1=1-22*x0-3*x0^2, p2=-11-3*x0, p3=-1); my(q0=1-22*x0-3*x0^2, q1=-22-6*x0, q2=-3); my(r0=-(3+x0), r1=-1); my(c=vector(MT+3)); c[1]=y0; c[2]=y1; for(n=0, MT, my(cn1=c[n+2], cn=c[n+1], cm1=if(n>=1, c[n], 0)); my(rhs = Rv[n+1] - (p1*(n+1)*n*cn1 + p2*n*(n-1)*cn + p3*(n-1)*(n-2)*cm1 + q0*(n+1)*cn1 + q1*n*cn + q2*(n-1)*cm1 + r0*cn + r1*cm1)); c[n+3] = rhs/(p0*(n+2)*(n+1))); my(s=0, sd=0); forstep(k=MT+2, 0, -1, s = s*h + c[k+1]); forstep(k=MT+2, 1, -1, sd = sd*h + k*c[k+1]); [s, sd];
rad(x) = my(d=[abs(x), abs(x-t1), abs(x-t2)]); vecmin(d);
/* adaptive walk from x0 to x1 along the straight segment */
walk(st, x0, x1, kind, par) = my(y=st[1], yp=st[2], cur=x0); while(abs(x1-cur) > 1e-60, my(rr=rad(cur)); my(rem=abs(x1-cur)); my(h=min(rem, 0.4*rr)); my(dir=(x1-cur)/rem); my(Rv = if(kind==0, rser0(cur,MT), if(kind==1, rserC(cur,MT,1), rserP(cur, MT, par[1], par[2])))); my(z=stp(y, yp, cur, h*dir, Rv)); y=z[1]; yp=z[2]; cur = cur + h*dir); [y, yp];
/* path = list of vertices; returns [y,y'] at the last vertex */
run(cf, kind, par, verts) = my(x0=verts[1]); my(st=[ev(cf,x0), evd(cf,x0)]); for(j=2, #verts, st = walk(st, verts[j-1], verts[j], kind, par)); st;
/* ---------- loop around t2 ---------- */
Nc = 32; rho = 5.0; xb2 = t2 + rho;
appr2 = [-0.05, -0.12, -0.25, -0.5, -0.9, -1.5, -2.4, -3.6, -5.0, xb2];
circ2 = vector(Nc+1, i, t2 + rho*exp(2*Pi*I*(i-1)/Nc));
path2 = concat(appr2, circ2[2..Nc+1]);
print("");
print("=== monodromy around t2 = ", t2, " (base point ", xb2, ") ===");
gA  = run(Av,   0, 0, path2);
gD  = run(BDv,  1, 0, path2);
gNp = run(Bnwp, 2, [2, ph5], path2);
gN  = run(Bnw,  2, [2, -phm5], path2);
sA = [ev(Av,xb2), evd(Av,xb2)]; sD = [ev(BDv,xb2), evd(BDv,xb2)];
/* the start values at xb2 must be obtained by continuation too */
st2A = run(Av, 0, 0, appr2); st2D = run(BDv, 1, 0, appr2); st2Np = run(Bnwp, 2, [2,ph5], appr2); st2N = run(Bnw, 2, [2,-phm5], appr2);
dA = gA[1]-st2A[1]; dD = gD[1]-st2D[1]; dNp = gNp[1]-st2Np[1]; dN = gN[1]-st2N[1];
dAp = gA[2]-st2A[2]; dDp = gD[2]-st2D[2]; dNpp = gNp[2]-st2Np[2];
print("  Delta A       = ", dA);
print("  Delta B_D     = ", dD);
print("  Delta B'_new  = ", dNp);
print("  Delta B_new   = ", dN, "   (expected: NOT a multiple of u -- log^2 at t2)");
print("  pi_D  = Delta B_D / Delta A     = ", dD/dA);
print("  pi'   = Delta B'_new / Delta A  = ", dNp/dA);
print("  consistency (derivatives): dD'/dA' = ", dDp/dAp, "   dNp'/dA' = ", dNpp/dAp);
write("far.txt", "piD = ", dD/dA);
write("far.txt", "pip = ", dNp/dA);
write("far.txt", "piDd = ", dDp/dAp);
write("far.txt", "pipd = ", dNpp/dAp);
quit;
