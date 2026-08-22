
/* 20_scan.c -- exhaustive search of Calegari's *shifted* order-2 normalisation
 *
 *    (n+1)^2 u_{n+1} = P(n) u_n - c (n-r)^2 u_{n-1},   P(n) = aq n^2 + bq n + cq
 *
 * for rows whose A-solution (A_{-1}=0, A_0=1, recurrence from n=0) is INTEGRAL.
 * Calegari's 2-adic Catalan row (math/0408214 sec.4) is (-32,0,4,256,1);
 * Zagier's normalisation is r=0, bq=aq (e.g. Zagier C = (10,10,3,9,0)).
 *
 * Only rows that could possibly score S_p > -0.3 are emitted:
 *      S_p = v_p(c) log p - k - log lambda_1,  lambda_1 = larger |char. root|,
 * so with k >= 1 we need lambda_1 <= exp(v_p(c) log p - 0.7), which bounds
 *      |aq| <= lambda_1 + |c|/lambda_1 .
 * (k = 0 -- an integral second solution -- is NOT covered; see the write-up.)
 *
 * Integrality is tested exactly in __int128 until the values overflow (>= 12
 * verified divisibility steps are required, which already leaves only structured
 * families); every survivor is re-verified exactly to n = 300 by 20_verify.gp.
 *
 * Build: gcc -O2 -fopenmp -o .../20_scan .../20_scan.c
 * Run:   .../20_scan > .../20_scan.out
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef __int128 i128;

static const long long CLIST[] = {
    64, 128, 256, 27, 81, 243, 25, 125, 49,
   -64,-128,-256,-27,-81,-243,-25,-125,-49
};
#define NC (int)(sizeof(CLIST)/sizeof(CLIST[0]))

#define BLIM 3000
#define CLIM 3000      /* cq >= 0 only: (aq,bq,cq) ~ (-aq,-bq,-cq) via u_n -> (-1)^n u_n */
#define RMAX 3
#define NA    40
#define NSTEP 12       /* required number of verified divisibility steps */
#define KMIN  1.0      /* smallest denominator exponent considered */
#define MARGIN 0.3     /* report threshold S_p > -MARGIN */

/* v_p(c) log p for the unique prime p dividing |c| */
static double vplogp(long long c)
{
    long long a = c < 0 ? -c : c;
    for (long long q = 2; q <= a; q++)
        if (a % q == 0) { int v = 0; while (a % q == 0) { a /= q; v++; } return v * log((double)q); }
    return 0.0;
}

int main(void)
{
    printf("# aq bq cq c r nsteps   [(n+1)^2 u_{n+1} = (aq n^2+bq n+cq) u_n - c (n-r)^2 u_{n-1}]\n");
    fflush(stdout);
    for (int ci = 0; ci < NC; ci++) {
        long long c = CLIST[ci];
        double lammax = exp(vplogp(c) - KMIN + MARGIN);
        long long ac = c < 0 ? -c : c;
        int alim = (int)floor(lammax + (double)ac / lammax) + 1;
        fprintf(stderr, "c=%lld  lambda_1 <= %.3f  |aq| <= %d\n", c, lammax, alim);
        for (int r = 0; r <= RMAX; r++) {
            #pragma omp parallel for schedule(dynamic,1)
            for (int aq = -alim; aq <= alim; aq++) {
                /* archimedean dominant root from x^2 - aq x + c */
                double disc = (double)aq * aq - 4.0 * (double)c;
                double lam1 = (disc >= 0) ? (fabs((double)aq) + sqrt(disc)) / 2.0 : sqrt((double)ac);
                if (lam1 > lammax) continue;
                for (long long bq = -BLIM; bq <= BLIM; bq++) {
                    for (long long cq = 0; cq <= CLIM; cq++) {
                        long long p1 = aq + bq + cq;
                        long long s1 = (long long)(1 - r) * (1 - r);
                        long long t1 = p1 * cq - c * s1;
                        if (t1 & 3) continue;
                        long long a2 = t1 >> 2;
                        long long p2 = 4LL * aq + 2LL * bq + cq;
                        long long s2 = (long long)(2 - r) * (2 - r);
                        long long t2 = p2 * a2 - c * s2 * cq;
                        if (t2 % 9) continue;
                        long long a3 = t2 / 9;
                        i128 prev = a2, cur = a3;
                        int n, ok = 1;
                        const i128 CAP = ((i128)1) << 100;
                        for (n = 3; n < NA; n++) {
                            i128 pn = (i128)aq * n * n + (i128)bq * n + cq;
                            i128 sh = (i128)(n - r) * (n - r);
                            i128 t  = pn * cur - (i128)c * sh * prev;
                            i128 d  = (i128)(n + 1) * (n + 1);
                            if (t % d) { ok = 0; break; }
                            prev = cur; cur = t / d;
                            if (cur > CAP || cur < -CAP) break;
                        }
                        if (!ok || n < NSTEP) continue;
                        if (cur == 0 || prev == 0) continue;   /* degenerate / terminating */
                        #pragma omp critical
                        { printf("%d %lld %lld %lld %d %d\n", aq, bq, cq, c, r, n); }
                    }
                }
            }
            fflush(stdout);
        }
    }
    return 0;
}
