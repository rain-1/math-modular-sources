default(parisizemax, 8000000000);
default(realprecision, 1500);
read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/17_arith_lib.gp");
fn = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/arith_arch.log";
logit(s) = { print(s); write(fn, s); };
NN = 300;
ns = [50,100,150,200,250,300];
{DN = vector(NN+1); DN[1]=1; for(n=1,NN, DN[n+1]=lcm(DN[n],n));}
z2 = Pi^2/6;
lm3 = lfun(-3,2);
{XI = [Catalan/4, Catalan/2 - 3*z2/16, 3*z2/8 - Catalan/2,
       (2*z2 + 15*lm3)/32, (15*lm3 - 6*z2)/16, (2*z2 - 3*lm3)/8];}
{XN = ["G/4","G/2-3z(2)/16","3z(2)/8-G/2","(2z(2)+15L(2,chi-3))/32",
       "(15L(2,chi-3)-6z(2))/16","(2z(2)-3L(2,chi-3))/8"];}

arch(nm, xi, xin, av, bv) = {
  my(ff, rr, ra, chk);
  logit(Str("--- ", nm, "   xi = ", xin, " = ", precision(xi,30), " ---"));
  chk = precision(bv[NN+1]/av[NN+1] - xi, 20);
  logit(Str("  b_N/a_N - xi at N=", NN, " : ", chk));
  logit(Str("  n : (1/n)log|d_n^2(xi a_n - b_n)|   (1/n)log a_n   (1/n)log|xi a_n - b_n|"));
  for(j=1,#ns, my(n=ns[j], lf = xi*av[n+1]-bv[n+1]);
    ff = DN[n+1]^2*lf;
    logit(Str("  ", n, " : ", precision(log(abs(ff))/n,12), "   ",
              precision(log(abs(1.0*av[n+1]))/n,12), "   ",
              precision(log(abs(lf))/n,12))));
  logit("");
};

{for(i=1,#ROWS, my(rr = row4(ROWS[i][2], NN)); arch(ROWS[i][1], XI[i], XN[i], rr[1], rr[2]));}
{my(rr = zrow2(12,4,32,NN)); arch("Zagier E (12,4,32)", Catalan/2, "G/2", rr[1], rr[2]);}
logit("DONE");
quit;
