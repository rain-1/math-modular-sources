default(parisizemax, 8000000000);
default(threadsizemax, 4000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/census_util.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/multi_prime/out/census_2_nest.log";
NMAX = 60;

/* Nesterenko (4,7) Catalan row, exact partial fractions.
   R_n(t) = c * prod_{k=0}^{4n-1}(t-k) / prod_{r=0}^{3n}(t+r+1/2)^2,
   c = (3n)! prod_{i=0}^{3n-2}(3/2+i) / (4n)!.
   raw linear form:  J_n = 4 B_n * Catalan - C_n,  so Qraw = 4 B_n, Praw = C_n.  */
{
nestpf(n) =
 my(NB=3*n, cc=(3*n)!*prod(i=0,3*n-2,3/2+i)/(4*n)!,
    num=cc*prod(k=0,4*n-1,('t-k)), tail=vector(NB+2), pre=vector(NB+2), mat, rhs, sol);
 pre[1]=1; for(j=1,NB, pre[j+1]=pre[j]*('t+j));
 tail[NB+2]=1; forstep(j=NB,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 my(bs=vector(2*(NB+1)), lin=vector(NB+2));
 lin[1]=1; for(j=0,NB, lin[j+2]=lin[j+1]*('t+j+1/2));
 for(j=0,NB, bs[2*j+1]=pre[j+1]*lin[j+2]*tail[j+2];
             bs[2*j+2]=pre[j+1]*lin[j+1]*tail[j+2]);
 my(deg=6*n+1, KK=2*(NB+1));
 mat = matrix(deg+1, KK, r, c2, polcoeff(bs[c2], r-1, 't));
 rhs = vectorv(deg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(mat, rhs);
 [vector(NB+1,j,sol[2*(j-1)+1]), vector(NB+1,j,sol[2*(j-1)+2])];
}
{
nestBC(n) = my(NB=3*n, r=nestpf(n), a1=r[1], a2=r[2], bb, cc2, gg);
 bb = sum(j=0,NB, a2[j+1]);
 gg = vector(NB+1); gg[1]=1;
 for(v=1,NB, gg[v+1] = gg[v]*v/(3/2+v-1));
 cc2 = sum(j=1,NB, sum(v=0,j-1, gg[v+1]*(a1[j+1]+a2[j+1]/(v+1/2))));
 [bb, cc2];
}

qv = vector(NMAX+1); pv = vector(NMAX+1);
{ my(t0=getabstime());
  for(n=1,NMAX, my(r=nestBC(n)); qv[n+1] = 4*r[1]; pv[n+1] = r[2];
    if(n % 10 == 0, print("  built n=",n,"  ms=",getabstime()-t0)));
  qv[1] = 4*qv[2]/qv[2]*0 + 1; pv[1] = 0; }
/* index 0 slot is a dummy; all reporting uses n>=2 */

logit(fn, Str("### Nesterenko (4,7) Catalan row, n=1..", NMAX, " ###"));
logit(fn, Str("  Qraw_n = 4 B_n, Praw_n = C_n, with J_n = Qraw_n*Catalan - Praw_n"));
/* cross-check against lattice/catalan_audit/nest_exact.txt  (V_n = 4^{7n+1} D^2 B_n, U_n = 4^{7n} D^2 C_n) */
{ my(ok=1, dd, vv, uu, ex=[[49575859200,45409750528],
      [183916572822218396160000,168461252904103170816000],
      [827587030023457009359480573898752000,758041245688795791789191614048845824],
      [4088867766240839623814891272419112710329026560000,3745262193016868931308061746500004874112221634560]]);
  for(n=1,4, dd = lcm(vector(6*n,i,i));
    vv = 4^(7*n)*dd^2*qv[n+1]; uu = 4^(7*n)*dd^2*pv[n+1];
    if(vv != ex[n][1] || uu != ex[n][2], ok=0));
  logit(fn, Str("  cross-check vs catalan_audit/nest_exact.txt (V_n,U_n), n=1..4: ", ok)); }
{ logit(fn, Str("  numeric U/(V*Catalan) at n=8 : ", 1.0*pv[9]/(qv[9]*Catalan))); }

/* denominator structure of Qraw and Praw */
logit(fn, "--- denominator factorisations ---");
{ my(dq, dp);
  for(n=1,NMAX, if(n<=6 || n==20 || n==40 || n==NMAX,
    dq = denominator(qv[n+1]); dp = denominator(pv[n+1]);
    logit(fn, Str("  n=",n,"  den(Qraw)=2^",valuation(dq,2)," oddpart=",dq/2^valuation(dq,2),
                  "   v2(den Praw)=",valuation(dp,2)," oddpart(den Praw)=",dp/2^valuation(dp,2))))); }

logit(fn, "--- kappa: v_p(den Qraw_n) for all p<=50 ---");
{ my(pl = primes(15), s);
  for(i=1,#pl, s = "";
    for(k=1,5, my(n=[10,20,30,45,NMAX][k]); s = Str(s," n=",n,":",valuation(denominator(qv[n+1]),pl[i])));
    logit(fn, Str("  p=",pl[i],":",s))); }

ns = [15, 30, 45, NMAX-1, NMAX];
report(fn, "Nesterenko(4,7)", qv, pv, primes(15), ns);

logit(fn, "--- fine v_2 data, n=NMAX-15..NMAX ---");
{ my(s1="", s2="", s3="");
  for(n=NMAX-15,NMAX,
    s1 = Str(s1," ",valuation(denominator(qv[n+1]),2));
    if(n>1, s2 = Str(s2," ",valuation(pv[n+1]/qv[n+1]-pv[n]/qv[n],2))));
  logit(fn, Str("  v2(den Qraw):",s1));
  logit(fn, Str("  v2(incr)    :",s2)); }
logit(fn, "DONE");
quit;
