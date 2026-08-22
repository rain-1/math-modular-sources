/* Beukers-type square-root rows: exhaustive integrality scan of the
 * R3-normalised second-order class
 *
 *   (n+1)^2 a_{n+1} = (al n^2 + be n + ga) a_n - (de n^2 + ep n + ze) a_{n-1},
 *   a_0 = 1, a_1 = ga.
 *
 * The R3 (free-integration) normalisation Q = (1/2) theta P forces
 *   be = al/2,  ep = -de,   ze free
 * and gives local exponents (0,1/2) at both finite singular points.
 * Sub-case ze = de/4 is the "sqrt class" of SQRT_APERY Thm 2.
 * Zagier class (exponents (0,0)) is be = al, ep = 0, ze = 0.
 *
 * Characteristic roots: x^2 - al x + de = 0.   score = log(1/|lam2|) - 2.
 *
 * Integrality is tested division-free via A_n = a_n (n!)^2:
 *   A_{n+1} = P(n) A_n - n^2 Q(n) A_{n-1},   a_n in Z  <=>  v_p(A_n) >= 2 v_p(n!).
 * Each prime is tested separately modulo p^K with K chosen so no precision is lost.
 *
 * Usage: 01_class_scan  CLASS ALMIN ALMAX DEMAX GAMAX ZEMAX NTEST
 *   CLASS = 0 : R3 class (be=al/2, ep=-de, ze free)
 *   CLASS = 1 : sqrt class (ze = de/4)  [3 parameters]
 *   CLASS = 2 : Zagier class (be=al, ep=0, ze=0) [3 parameters]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef __int128 i128;
typedef unsigned __int128 u128;

static int NP, PR[32], KK[32];
static u128 MOD[32];
static int V2FACT[256], VPF[32][256];

/* v_p(n!) */
static int vpfact(int n,int p){int s=0;long long q=p;while(q<=n){s+=n/q; if(q> (long long)n/p) break; q*=p;} return s;}

/* run the A-recursion mod p^K, return 1 if v_p(A_n) >= 2 v_p(n!) for all n<=N */
static inline int testp(long long al,long long be,long long ga,
                        long long de,long long ep,long long ze,
                        int pi,int N)
{
    u128 M = MOD[pi];
    int p = PR[pi];
    u128 A0=1, A1;
    long long g = ga;
    /* reduce to [0,M) */
    #define RED(x) ( (x)>=0 ? (u128)(x)%M : M - ((u128)(-(x))%M) )
    A1 = RED(g);
    /* n=1 check trivially ok (v_p(1!)=0) */
    for(int n=1;n<=N-1;n++){
        long long Pn = al*(long long)n*n + be*(long long)n + ga;
        long long Qn = de*(long long)n*n + ep*(long long)n + ze;
        u128 pP = RED(Pn), pQ = RED(Qn);
        u128 nn = (u128)((long long)n*n) % M;
        u128 t1 = (pP*A1) % M;
        u128 t2 = (((nn*pQ)%M)*A0) % M;
        u128 A2 = (t1 + M - t2) % M;
        /* need v_p(A2) >= 2 v_p((n+1)!) */
        int need = 2*VPF[pi][n+1];
        if(need>0){
            if(A2==0){ /* fine, but keep going with zero */ }
            else{
                u128 x=A2; int v=0;
                while(v<need && x % (u128)p == 0){ x /= (u128)p; v++; }
                if(v<need) return 0;
            }
        }
        A0=A1; A1=A2;
    }
    return 1;
    #undef RED
}

int main(int argc,char**argv)
{
    int CLASS = atoi(argv[1]);
    long long ALMIN=atoll(argv[2]), ALMAX=atoll(argv[3]);
    long long DEMAX=atoll(argv[4]), GAMAX=atoll(argv[5]), ZEMAX=atoll(argv[6]);
    int N = atoi(argv[7]);

    /* primes and moduli */
    int prs[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,0};
    NP=0;
    for(int i=0;prs[i];i++){
        int p=prs[i]; if(p>N) break;
        int need = 2*vpfact(N,p) + 2;
        /* modulus p^K <= 2^126 */
        u128 M=1; int K=0;
        while(K<need){ if(M > (((u128)1)<<126)/p) break; M*= (u128)p; K++; }
        if(K<need) { fprintf(stderr,"precision fail p=%d\n",p); }
        PR[NP]=p; KK[NP]=K; MOD[NP]=M;
        for(int n=0;n<=N+1;n++) VPF[NP][n]=vpfact(n,p);
        NP++;
    }
    fprintf(stderr,"# primes tested: %d (up to %d), N=%d\n",NP,PR[NP-1],N);

    long long count=0, hits=0;
    for(long long al=ALMIN; al<=ALMAX; al++){
        if(CLASS!=2 && (al%2)) continue;
        for(long long de=-DEMAX; de<=DEMAX; de++){
            if(de==0) continue;
            long long be = (CLASS==2)? al : al/2;
            long long ep = (CLASS==2)? 0  : -de;
            long long zemin, zemax, zestep;
            if(CLASS==1){ if(de%4) continue; zemin=de/4; zemax=de/4; zestep=1; }
            else if(CLASS==2){ zemin=0; zemax=0; zestep=1; }
            else { zemin=-ZEMAX; zemax=ZEMAX; zestep=1; }
            for(long long ga=-GAMAX; ga<=GAMAX; ga++){
                for(long long ze=zemin; ze<=zemax; ze+=zestep){
                    /* fast n=1 gate: A_2 = P(1)*ga - Q(1); need 4 | A_2 */
                    long long A2 = (al+be+ga)*ga - (de+ep+ze);
                    if(A2 & 3) continue;
                    count++;
                    int ok=1;
                    for(int i=0;i<NP;i++){ if(!testp(al,be,ga,de,ep,ze,i,N)){ok=0;break;} }
                    if(ok){ hits++; printf("%lld %lld %lld %lld %lld %lld\n",al,be,ga,de,ep,ze); }
                }
            }
        }
    }
    fprintf(stderr,"# gate-passed %lld, hits %lld\n",count,hits);
    return 0;
}
