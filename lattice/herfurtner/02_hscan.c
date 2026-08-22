/* Herfurtner-class integrality scan.
 *
 * Row:  (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1},  u_0 = 1, u_1 = P(0).
 *   P(n) = A n^2 + A(1-rho) n + B            (rho = (j1+j2)/(2M))
 *   Q(n) = C (M n - j1)(M n - j2)
 *          = C M^2 n^2 - C M (j1+j2) n + C j1 j2
 *
 * The class (M; j1, j2) fixes the local exponents:
 *   exponents at the two finite singular points : (0, rho), (0, rho)
 *   exponents at infinity                       : (1 - j1/M, 1 - j2/M)
 *   delta_infty = (j2 - j1)/M
 * A(1-rho) must be an integer -> A must be divisible by 2M/gcd(2M, 2M-j1-j2).
 *
 * Integrality is tested division-free on U_n = u_n (n!)^2:
 *   U_{n+1} = P(n) U_n - n^2 Q(n) U_{n-1},   u_n in Z  <=>  v_p(U_n) >= 2 v_p(n!).
 * Each prime p <= N is tested separately modulo p^K.
 *
 * usage: 02_hscan M J1 J2 AMIN AMAX BMAX CMAX N
 * output lines: M J1 J2 A B C
 */
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
static int NP, PR[32];
static u128 MOD[32];
static int VPF[32][300];

static int vpfact(int n,int p){int s=0;long long q=p;while(q<=n){s+=n/q;if(q>(long long)n/p)break;q*=p;}return s;}

static inline int testp(long long al,long long be,long long ga,
                        long long de,long long ep,long long ze,int pi,int N)
{
    u128 M = MOD[pi]; int p = PR[pi];
    u128 A0=1, A1;
    #define RED(x) ( (x)>=0 ? (u128)(x)%M : M - ((u128)(-(x))%M) )
    A1 = RED(ga);
    for(int n=1;n<=N-1;n++){
        long long Pn = al*(long long)n*n + be*(long long)n + ga;
        long long Qn = de*(long long)n*n + ep*(long long)n + ze;
        u128 pP = RED(Pn), pQ = RED(Qn);
        u128 nn = (u128)((long long)n*n) % M;
        u128 A2 = ((pP*A1)%M + M - (((nn*pQ)%M)*A0)%M) % M;
        int need = 2*VPF[pi][n+1];
        if(need>0 && A2){
            u128 x=A2; int v=0;
            while(v<need && x % (u128)p == 0){ x /= (u128)p; v++; }
            if(v<need) return 0;
        }
        A0=A1; A1=A2;
    }
    return 1;
    #undef RED
}

int main(int argc,char**argv)
{
    long long M=atoll(argv[1]), J1=atoll(argv[2]), J2=atoll(argv[3]);
    long long AMIN=atoll(argv[4]), AMAX=atoll(argv[5]);
    long long BMAX=atoll(argv[6]), CMAX=atoll(argv[7]);
    int N=atoi(argv[8]);

    int prs[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,0};
    NP=0;
    for(int i=0;prs[i];i++){
        int p=prs[i]; if(p>N) break;
        int need = 2*vpfact(N,p)+2;
        u128 Mo=1; int K=0;
        while(K<need){ if(Mo > (((u128)1)<<126)/p) break; Mo*=(u128)p; K++; }
        if(K<need) fprintf(stderr,"# precision short p=%d (%d<%d)\n",p,K,need);
        PR[NP]=p; MOD[NP]=Mo;
        for(int n=0;n<=N+1;n++) VPF[NP][n]=vpfact(n,p);
        NP++;
    }
    /* A must satisfy A*(2M-J1-J2) % (2M) == 0 */
    long long S = 2*M-J1-J2;
    long long g = 2*M, x = S<0?-S:S, y=g;
    while(x){ long long r=y%x; y=x; x=r; }   /* y = gcd(|S|,2M) */
    long long Astep = (2*M)/y;
    long long a0 = AMIN;
    if(a0 % Astep) a0 += Astep - ((a0 % Astep)+Astep)%Astep;
    fprintf(stderr,"# class M=%lld j=(%lld,%lld) Astep=%lld primes=%d N=%d\n",M,J1,J2,Astep,NP,N);

    long long hits=0;
    for(long long A=a0; A<=AMAX; A+=Astep){
        long long be = A*S/(2*M);
        for(long long C=-CMAX; C<=CMAX; C++){
            if(C==0) continue;
            long long de = C*M*M, ep = -C*M*(J1+J2), ze = C*J1*J2;
            for(long long B=-BMAX; B<=BMAX; B++){
                if(B==0) continue;                      /* u_1 = 0 -> degenerate */
                long long U2 = (A+be+B)*B - (de+ep+ze); /* need 4 | U2 */
                if(U2 & 3) continue;
                int ok=1;
                for(int i=0;i<NP;i++) if(!testp(A,be,B,de,ep,ze,i,N)){ok=0;break;}
                if(ok){ hits++; printf("%lld %lld %lld %lld %lld %lld\n",M,J1,J2,A,B,C); }
            }
        }
    }
    fprintf(stderr,"# hits %lld\n",hits);
    return 0;
}
