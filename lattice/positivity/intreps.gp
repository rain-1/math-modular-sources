/* lattice/positivity/intreps.gp  --- Task 1: independent numerical verification
   of the positive-kernel integral representations, at >= 30 digits.
   Standalone (does not need rows_pos.gp except for `mom`).                  */

/* Catalan-world Beukers moment, evaluated as an all-positive series:
   momnum(m,j) = sum_k binom(k+a,a) B(a+1/2+k,b+1) B(a+1+k,b+1/2),  a=m+j,b=m
   = int_0^1 int_0^1 x^{a-1/2}(1-x)^b y^a (1-y)^{b-1/2} (1-xy)^{-(a+1)}.
   Every term is positive, so the value is manifestly > 0.                    */
{
momnum(m,jj) = my(a=m+jj, b=m, tk, tot, eps, k);
 tk = gamma(a+1/2)*gamma(b+1)/gamma(a+b+3/2) * gamma(a+1)*gamma(b+1/2)/gamma(a+b+3/2);
 tot = tk; k = 0; eps = 10.^(-(precision(1.)-8));
 while(tk > tot*eps && k < 300000,
   tk *= (k+1+a)/(k+1) * (a+1/2+k)/(a+b+3/2+k) * (a+1+k)/(a+b+3/2+k);
   tot += tk; k++);
 [tot, k];
}

/* Beukers' zeta(2) kernel, index n:  all-positive series                    */
/* Beukers' zeta(2) kernel int int (x(1-x)y(1-y))^n/(1-xy)^{n+1}: all-positive */
{
beuk2(n) = my(tk, tot, eps, k);
 tk = (gamma(n+1)^2/gamma(2*n+2))^2;
 tot = tk; k = 0; eps = 10.^(-(precision(1.)-8));
 while(tk > tot*eps && k < 300000,
   tk *= (k+1+n)/(k+1) * ((n+1+k)/(2*n+2+k))^2;
   tot += tk; k++);
 [tot, k];
}

/* Beukers' zeta(3) triple integral, index n, as an iterated numeric integral */
{
beuk3(n) = intnum(x=0,1, intnum(y=0,1, intnum(z=0,1,
   (x*(1-x)*y*(1-y)*z*(1-z))^n / (1-(1-x*y)*z)^(n+1))));
}

/* L(2,chi_{-3}) Euler kernel  (ONE_CLASS_TWO_WORLDS Thm 1 remark)           */
chi3ker() = intnum(x=0,1, intnum(y=0,1, 1/(1+x*y+x^2*y^2)));
