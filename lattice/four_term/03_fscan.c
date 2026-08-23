/* Four-term (rank 2, five singular point) integrality scan.
 *
 * Row:  (n+1)^2 u_{n+1} = P(n)u_n - Q(n)u_{n-1} + R(n)u_{n-2},  u_0=1, u_{-1}=u_{-2}=0
 *   P(n) = a n^2 + a(1-rho) n + c
 *   Q(n) = d n^2 - 2 rho d n + f
 *   R(n) = C (M n - J1)(M n - J2)      with  J1+J2 = (1+3rho) M
 * Class (rho; M, J1, J2) fixes the local exponents:
 *   (0,0) at t=0 ; (0,rho) at each root t_i of  1 - a t + d t^2 - g t^3 ;
 *   (2-J1/M, 2-J2/M) at infinity.
 *
 * Integrality tested division-free on U_n = u_n (n!)^2:
 *   U_{n+1} = P(n) U_n - n^2 Q(n) U_{n-1} + n^2 (n-1)^2 R(n) U_{n-2},
 *   u_n in Z  <=>  v_p(U_n) >= 2 v_p(n!)     for every prime p.
 *
 * usage: 03_fscan RN RD M J1 J2 AMIN AMAX CMAX DMAX FMAX GMAX N [ASTRIDE ASHIFT]
 *        rho = RN/RD.
 * out:   RN RD M J1 J2 a c d f C
 */
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
typedef long long ll;
static int NP, PR[32];
static u128 MOD[32];
static int VPF[32][320];

#define N2 14
#define K2 24
#define MSK2 ((1ULL<<K2)-1)
static int VP2[N2+2];
static unsigned long long WM2[N2+2], C1M[N2+2], FCM[N2+2];
/* mask-arithmetic 2-adic prefilter; P,Q,W precomputed, only C varies */
static inline int test2(ll al,ll be,ll ga,ll de,ll ep,ll ze,unsigned long long Cm)
{
    unsigned long long r0=0, r1=1, r2=((unsigned long long)ga)&MSK2;
    for(int n=1;n<N2;n++){
        ll nn=n;
        unsigned long long pP=((unsigned long long)(al*nn*nn+be*nn+ga))&MSK2;
        unsigned long long pQ=((unsigned long long)(de*nn*nn+ep*nn+ze))&MSK2;
        unsigned long long pR=(Cm*WM2[n])&MSK2;
        unsigned long long nx = (pP*r2)&MSK2;
        nx = (nx - ((C1M[n]*pQ)&MSK2)*r1)&MSK2;
        nx = (nx + ((FCM[n]*pR)&MSK2)*r0)&MSK2;
        int need = 2*VP2[n+1];
        if(need && nx && (int)__builtin_ctzll(nx) < need) return 0;
        r0=r1;r1=r2;r2=nx;
    }
    return 1;
}

static int NS=12, NPS, PRS[8], KS[8];
static unsigned long long MODS[8];
static int VPFS[8][20];

/* fast 64-bit prefilter: primes 2,3,5,7,11 to n = NS, modulus p^K < 2^31 */
static inline int test_small(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io)
{
    for(int i=0;i<NPS;i++){
        unsigned long long MM = MODS[i]; int p = PRS[i];
        unsigned long long r0=0, r1=1, r2;
        { ll v = ga % (ll)MM; r2 = (unsigned long long)(v<0? v+(ll)MM : v); }
        for(int n=1;n<=NS-1;n++){
            ll nn=n;
            ll Pn = (al*nn*nn + be*nn + ga) % (ll)MM;
            ll Qn = (de*nn*nn + ep*nn + ze) % (ll)MM;
            ll Rn = (et*nn*nn + th*nn + io) % (ll)MM;
            unsigned long long pP=(unsigned long long)(Pn<0?Pn+(ll)MM:Pn);
            unsigned long long pQ=(unsigned long long)(Qn<0?Qn+(ll)MM:Qn);
            unsigned long long pR=(unsigned long long)(Rn<0?Rn+(ll)MM:Rn);
            unsigned long long c1=(unsigned long long)(nn*nn)%MM;
            unsigned long long fc=(unsigned long long)(nn*(nn-1))%MM; fc=(fc*fc)%MM;
            unsigned long long nx = (pP*r2)%MM;
            nx = (nx + MM - (((c1*pQ)%MM)*r1)%MM)%MM;
            nx = (nx + (((fc*pR)%MM)*r0)%MM)%MM;
            int need = 2*VPFS[i][n+1];
            if(need>0 && nx){
                unsigned long long x=nx; int v=0;
                while(v<need && x%(unsigned long long)p==0){x/=(unsigned long long)p;v++;}
                if(v<need) return 0;
            }
            r0=r1;r1=r2;r2=nx;
        }
    }
    return 1;
}

static int vpfact(int n,int p){int s=0;ll q=p;while(q<=n){s+=n/q;if(q>(ll)n/p)break;q*=p;}return s;}
static ll gcdll(ll x,ll y){if(x<0)x=-x;if(y<0)y=-y;while(x){ll r=y%x;y=x;x=r;}return y;}

/* P = al n^2 + be n + ga ; Q = de n^2 + ep n + ze ; R = et n^2 + th n + io */
static inline int testp(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io,int pi,int N)
{
    u128 MM = MOD[pi]; int p = PR[pi];
    #define RED(x) ( (x)>=0 ? (u128)(x)%MM : MM - ((u128)(-(x))%MM) )
    u128 a0 = 0, a1 = 1, a2 = RED(ga);        /* U_{-1}, U_0, U_1 */
    for(int n=1;n<=N-1;n++){
        ll nn = n;
        ll Pn = al*nn*nn + be*nn + ga;
        ll Qn = de*nn*nn + ep*nn + ze;
        ll Rn = et*nn*nn + th*nn + io;
        u128 pP=RED(Pn), pQ=RED(Qn), pR=RED(Rn);
        u128 c1  = (u128)(nn*nn) % MM;
        u128 fac = (u128)(nn*(nn-1)) % MM;  fac = (fac*fac) % MM;   /* n^2 (n-1)^2 */
        u128 nx = (pP*a2) % MM;
        nx = (nx + MM - ((((c1*pQ)%MM)*a1) % MM)) % MM;
        nx = (nx + ((((fac*pR)%MM)*a0) % MM)) % MM;
        int need = 2*VPF[pi][n+1];
        if(need>0 && nx){
            u128 x=nx; int v=0;
            while(v<need && x % (u128)p == 0){ x /= (u128)p; v++; }
            if(v<need) return 0;
        }
        a0=a1; a1=a2; a2=nx;
    }
    return 1;
    #undef RED
}

static ll modinv(ll a, ll m){
    if(m==1) return 0;
    ll g=m, xx=0, x1=1, a1=((a%m)+m)%m;
    while(a1){ ll q=g/a1, t=g-q*a1; g=a1; a1=t; t=xx-q*x1; xx=x1; x1=t; }
    return ((xx%m)+m)%m;
}
/* A*C = B (mod q)  ->  C = r (mod m);  0 if unsolvable */
static inline int lincong(ll A, ll B, ll q, ll *r, ll *m){
    A=((A%q)+q)%q; B=((B%q)+q)%q;
    ll g=gcdll(A,q); if(B%g) return 0;
    ll q1=q/g; if(q1==1){*r=0;*m=1;return 1;}
    ll inv=modinv((A/g)%q1, q1);
    *r=(ll)(((__int128)((B/g)%q1)*inv)%q1); *m=q1; return 1;
}
static inline int crt2(ll r1,ll m1,ll r2,ll m2,ll *r,ll *m){
    ll g=gcdll(m1,m2); if(((r2-r1)%g+g)%g) return 0;
    ll m2g=m2/g, l=m1*m2g;
    ll inv=modinv((m1/g)%m2g, m2g);
    ll k=(ll)(((__int128)((((r2-r1)/g)%m2g+m2g)%m2g)*inv)%m2g);
    *r=(ll)((((__int128)r1+(__int128)m1*k)%l+l)%l); *m=l; return 1;
}

/* reject sequences that vanish identically from some point on */
static inline int nontrivial(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io)
{
    __int128 v[10]; v[0]=0; v[1]=0; v[2]=1;          /* U_{-2},U_{-1},U_0 */
    int z=0;
    for(int n=0;n<=6;n++){
        __int128 nx = (__int128)(al*(ll)n*n+be*n+ga)*v[n+2]
                    - (__int128)((ll)n*n)*(de*(ll)n*n+ep*n+ze)*v[n+1]
                    + (__int128)((ll)n*n*(n-1)*(n-1))*(et*(ll)n*n+th*n+io)*v[n];
        v[n+3]=nx;
        if(nx==0){ if(++z>=3) return 0; } else z=0;
    }
    return 1;
}

int main(int argc,char**argv)
{
    if(argc<13){fprintf(stderr,"usage: %s RN RD M J1 J2 AMIN AMAX CMAX DMAX FMAX GMAX N [STRIDE SHIFT]\n",argv[0]);return 1;}
    ll RN=atoll(argv[1]), RD=atoll(argv[2]);
    ll M=atoll(argv[3]), J1=atoll(argv[4]), J2=atoll(argv[5]);
    ll AMIN=atoll(argv[6]), AMAX=atoll(argv[7]);
    ll CMAX=atoll(argv[8]), DMAX=atoll(argv[9]), FMAX=atoll(argv[10]), GMAX=atoll(argv[11]);
    ll DMAX0=DMAX, GMAX0=GMAX;
    int N=atoi(argv[12]);
    ll STRIDE = (argc>13)?atoll(argv[13]):1, SHIFT=(argc>14)?atoll(argv[14]):0;

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
    for(int n=0;n<=N2+1;n++){ VP2[n]=vpfact(n,2);
        C1M[n]=((unsigned long long)((ll)n*n))&MSK2;
        unsigned long long fz=((unsigned long long)((ll)n*(n-1)))&MSK2; FCM[n]=(fz*fz)&MSK2; }
    NPS=0;
    for(int i=0;prs[i];i++){
        int p=prs[i]; if(p>NS) break;
        int need = 2*vpfact(NS,p)+2;
        unsigned long long Mo=1; int K=0;
        while(K<need && Mo <= (1ULL<<31)/p){ Mo*=(unsigned long long)p; K++; }
        if(K<need) continue;             /* skip primes needing too much precision */
        PRS[NPS]=p; MODS[NPS]=Mo; KS[NPS]=K;
        for(int n=0;n<=NS+1;n++) VPFS[NPS][n]=vpfact(n,p);
        NPS++;
    }
    fprintf(stderr,"# prefilter primes=%d up to n=%d\n",NPS,NS);
    /* b = a(1-rho) = a(RD-RN)/RD  integral  ->  a step */
    ll Astep = RD/gcdll(RD, RD-RN);
    /* e = -2 rho d = -2 RN d / RD  integral  ->  d step */
    ll Dstep = RD/gcdll(RD, 2*RN);
    /* check J1+J2 = (1+3rho)M   i.e.  RD*(J1+J2) = (RD+3RN)*M */
    if(RD*(J1+J2) != (RD+3*RN)*M){fprintf(stderr,"# BAD CLASS: J1+J2 != (1+3rho)M\n");return 2;}

    ll a0 = AMIN; if(a0%Astep) a0 += Astep - ((a0%Astep)+Astep)%Astep;
    fprintf(stderr,"# class rho=%lld/%lld M=%lld J=(%lld,%lld) Astep=%lld Dstep=%lld primes=%d N=%d stride=%lld/%lld\n",
            RN,RD,M,J1,J2,Astep,Dstep,NP,N,SHIFT,STRIDE);

    ll hits=0, cnt=0;
    ll W2 = (2*M-J1)*(2*M-J2);           /* R(2) = C*W2 */
    ll W3 = (3*M-J1)*(3*M-J2);           /* R(3) = C*W3 */
    for(int n=0;n<=N2+1;n++) WM2[n]=((unsigned long long)((M*(ll)n-J1)*(M*(ll)n-J2)))&MSK2;
    for(ll a=a0; a<=AMAX; a+=Astep){
        if(((a-a0)/Astep) % STRIDE != SHIFT) continue;
        ll be = a*(RD-RN)/RD;
        /* window modes:  -1 : |lambda_2| < 1        (|d|<=2a+5,   |g|<=a+2)
                         -2 : |lambda_2| < 1/e      (|d|<=0.74(a+1)+2, |g|<=0.136(a+1)+1) */
        ll DM, GM;
        if(DMAX0 == -1){ DM = 2*a+5; }
        else if(DMAX0 == -2){ DM = (74*(a+1))/100 + 2; }
        else DM = DMAX0;
        if(GMAX0 == -1){ GM = (a+2)/(M*M)+1; }
        else if(GMAX0 == -2){ GM = ((136*(a+1))/1000)/(M*M) + 1; }
        else GM = GMAX0;
        DMAX = DM; GMAX = GM;
        for(ll c=-CMAX; c<=CMAX; c++){
            ll P1 = a+be+c, P2 = 4*a+2*be+c, P3 = 9*a+3*be+c;
            ll d0=-DMAX; if(d0%Dstep) d0 += Dstep - ((d0%Dstep)+Dstep)%Dstep;
            for(ll d=d0; d<=DMAX; d+=Dstep){
                ll ep = -2*RN*d/RD;
                ll Q1a = d+ep, Q2a = 4*d+2*ep, Q3a = 9*d+3*ep;
                for(ll f=-FMAX; f<=FMAX; f++){
                    ll U2 = P1*c - (Q1a+f);            /* U_2 = P(1)U_1 - Q(1)U_0 */
                    if(U2 & 3) continue;               /* 4 | U_2                 */
                    ll A3 = 4*W2, K3 = P2*U2 - 4*(Q2a+f)*c;   /* U_3 = K3 + A3*C  */
                    ll r1,m1,r2,m2,rr,mm;
                    if(!lincong(A3, -K3, 36, &r1, &m1)) continue;      /* 36 | U_3  */
                    ll K4 = P3*K3 - 9*(Q3a+f)*U2;      /* U_4 = K4 + A4*C          */
                    ll A4 = P3*A3 + 36*c*W3;
                    if(!lincong(A4, -K4, 576, &r2, &m2)) continue;     /* 576 | U_4 */
                    if(!crt2(r1,m1,r2,m2,&rr,&mm)) continue;
                    ll Cst = -GMAX;
                    { ll t=((Cst-rr)%mm+mm)%mm; Cst-=t; if(Cst<-GMAX) Cst+=mm; }
                    for(ll C=Cst; C<=GMAX; C+=mm){
                        if(C==0) continue;
                        /* u_1=u_2=u_3=0  =>  u_n = 0 for all n >= 1 : trivial row */
                        if(c==0 && U2==0 && K3 + 4*W2*C == 0) continue;
                        /* P,Q,R all proportional to (n+1)^2 : constant-coefficient row */
                        if(a==c && be==2*a && d==f && ep==2*d
                           && C*J1*J2==C*M*M && C*M*(J1+J2)==-2*C*M*M) continue;
                        if(!test2(a,be,c,d,ep,f,((unsigned long long)C)&MSK2)) continue;
                        ll et=C*M*M, th=-C*M*(J1+J2), io=C*J1*J2;
                        {   /* discriminant of  x^3 - a x^2 + d x - g,  g = C M^2 */
                            __int128 D = (__int128)18*a*d*et - (__int128)4*a*a*a*et
                                       + (__int128)a*a*d*d - (__int128)4*d*d*d - (__int128)27*et*et;
                            if(D==0) continue;         /* repeated characteristic root */
                        }
                        if(!nontrivial(a,be,c,d,ep,f,et,th,io)) continue;
                        if(!test_small(a,be,c,d,ep,f,et,th,io)) continue;
                        cnt++;
                        int ok=1;
                        for(int i=0;i<NP;i++) if(!testp(a,be,c,d,ep,f,et,th,io,i,N)){ok=0;break;}
                        if(ok){ hits++; printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",RN,RD,M,J1,J2,a,c,d,f,C); }
                    }
                }
            }
        }
    }
    fprintf(stderr,"# deep %lld hits %lld\n",cnt,hits);
    return 0;
}
