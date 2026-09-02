{ligozat(n) = my(dv = divisors(n), k = #dv);
  matrix(k, k, i, j, my(c = dv[i], d = dv[j]);
    n/(24*gcd(c, n/c)*c) * gcd(c,d)^2/d );}
epsc(n,c) = eulerphi(gcd(c, n/c));
{show(n, r) = my(dv=divisors(n), A=ligozat(n), ord = A*r~);
  my(deg = sum(t=1,#dv, epsc(n,dv[t])*max(ord[t],0)));
  my(Cm2 = prod(t=1,#dv, (n/dv[t])^r[t]));
  print("N=",n," r=",r," div=",dv," ord=",Vec(ord~)," deg=",deg," C=",if(issquare(1/Cm2), sqrtint(1/Cm2), 1/Cm2));}
show(12, [-4,4,4,-4,-4,4]);        /* Domb */
show(6,  [-5,1,-1,5]);             /* Apery */
show(18, [-6,6,12,0,-12,-6,6,0,0]);  /* placeholder */
show(18, [-6,6,12,-12,-6,6]);      /* wrong length test */
quit;
