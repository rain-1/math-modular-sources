/* Four-term MIXED-exponent scan: rho_1 != rho_2 = rho_3.
 *
 * Theory (derived from FOUR_TERM_SCAN.md Thm F1, generalised):
 *   chi(lam) = lam^3 - a lam^2 + d lam - g = (lam - r)(lam^2 - s lam + p)
 *   W(lam)   = lam^3 + B1 lam^2 + B2 lam + B3 = chi(lam)*[1 + sum_i rho_i lam_i/(lam-lam_i)]
 *   (W is the reversal of T = Sc - t Rc'), i.e.
 *     W = chi + rho_p (lam - r)(s lam - 2 p) + rho_r r (lam^2 - s lam + p)
 *   giving, with a = r+s, d = s r + p, g = r p,
 *     b = (1-rho_r) r + (1-rho_p) s
 *     e = -rho_p (2 p + r s) - rho_r r s
 *     h = -(1 + 2 rho_p + rho_r) g
 *     j = C J1 J2,  C = g/M^2,  (J1+J2)/M = 1 + 2 rho_p + rho_r
 *   (the equal case rho_r = rho_p = rho reproduces b=(1-rho)a, e=-2 rho d,
 *    h=-(1+3rho)g.)
 *
 * A rational root r is forced: Cor F2.1 says rho is constant on Galois orbits,
 * so an unequal exponent needs a reducible cubic; and r in Q, chi monic in Z[lam]
 * => r in Z, and then s,p in Z by Gauss.
 *
 * usage: 03_fmix RPN RPD RRN RRD M J1 J2  RMIN RMAX AMAX CMAX DMAX FMAX GMAX N [STRIDE SHIFT]
 *   rho_p = RPN/RPD on the quadratic factor, rho_r = RRN/RRD on the rational root.
 * out:  RPN RPD RRN RRD M J1 J2 r a c d f C     (a = r+s, d = s r + p, g = C M^2 = r p)
 */
#include <stdio.h>
#include <stdlib.h>

typedef unsigned __int128 u128;
typedef long long ll;
static int NP, PR[32];
static u128 MOD[32];
static int VPF[32][320];

static int NS=12, NPS, PRS[8];
static unsigned long long MODS[8];
static int VPFS[8][20];

static int vpfact(int n,int p){int s=0;ll q=p;while(q<=n){s+=n/q;if(q>(ll)n/p)break;q*=p;}return s;}
static ll gcdll(ll x,ll y){if(x<0)x=-x;if(y<0)y=-y;while(x){ll r=y%x;y=x;x=r;}return y;}

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

static inline int testp(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io,int pi,int N)
{
    u128 MM = MOD[pi]; int p = PR[pi];
    #define RED(x) ( (x)>=0 ? (u128)(x)%MM : MM - ((u128)(-(x))%MM) )
    u128 a0 = 0, a1 = 1, a2 = RED(ga);
    for(int n=1;n<=N-1;n++){
        ll nn = n;
        ll Pn = al*nn*nn + be*nn + ga;
        ll Qn = de*nn*nn + ep*nn + ze;
        ll Rn = et*nn*nn + th*nn + io;
        u128 pP=RED(Pn), pQ=RED(Qn), pR=RED(Rn);
        u128 c1  = (u128)(nn*nn) % MM;
        u128 fac = (u128)(nn*(nn-1)) % MM;  fac = (fac*fac) % MM;
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

#define MP1 2305843009213693951ULL
#define MP2 2305843009213693739ULL
static inline unsigned long long mmul(unsigned long long x,unsigned long long y,unsigned long long m)
{ return (unsigned long long)(((u128)x*y)%m); }
static inline unsigned long long mred(ll x, unsigned long long m)
{ ll r = x % (ll)m; return (unsigned long long)(r<0? r+(ll)m : r); }
static inline int nontrivial(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io)
{
    unsigned long long MM[2]={MP1,MP2}; int zr[2]={0,0};
    for(int q=0;q<2;q++){
        unsigned long long m=MM[q], v0=0,v1=0,v2=1; int z=0;
        for(int n=0;n<=8;n++){
            unsigned long long pP=mred(al*(ll)n*n+be*n+ga,m);
            unsigned long long pQ=mred(de*(ll)n*n+ep*n+ze,m);
            unsigned long long pR=mred(et*(ll)n*n+th*n+io,m);
            unsigned long long c1=mred((ll)n*n,m);
            unsigned long long fc=mred((ll)n*(n-1),m); fc=mmul(fc,fc,m);
            unsigned long long nx=mmul(pP,v2,m);
            nx=(nx + m - mmul(mmul(c1,pQ,m),v1,m))%m;
            nx=(nx + mmul(mmul(fc,pR,m),v0,m))%m;
            if(nx==0){ if(++z>=3){ zr[q]=1; break; } } else z=0;
            v0=v1;v1=v2;v2=nx;
        }
    }
    return !(zr[0]&&zr[1]);
}

static ll modinv(ll a, ll m){
    if(m==1) return 0;
    ll g=m, xx=0, x1=1, a1=((a%m)+m)%m;
    while(a1){ ll q=g/a1, t=g-q*a1; g=a1; a1=t; t=xx-q*x1; xx=x1; x1=t; }
    return ((xx%m)+m)%m;
}
static inline int lincong(ll A, ll B, ll q, ll *r, ll *m){
    A=((A%q)+q)%q; B=((B%q)+q)%q;
    ll g=gcdll(A,q); if(g==0){ if(B) return 0; *r=0;*m=1;return 1; }
    if(B%g) return 0;
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

int main(int argc,char**argv)
{
    if(argc<16){fprintf(stderr,"usage: %s RPN RPD RRN RRD M J1 J2 RMIN RMAX AMAX CMAX DMAX FMAX GMAX N [STRIDE SHIFT]\n",argv[0]);return 1;}
    ll RPN=atoll(argv[1]), RPD=atoll(argv[2]), RRN=atoll(argv[3]), RRD=atoll(argv[4]);
    ll M=atoll(argv[5]), J1=atoll(argv[6]), J2=atoll(argv[7]);
    ll RMIN=atoll(argv[8]), RMAX=atoll(argv[9]), AMAX=atoll(argv[10]);
    ll CMAX=atoll(argv[11]), DMAX=atoll(argv[12]), FMAX=atoll(argv[13]), GMAX=atoll(argv[14]);
    int N=atoi(argv[15]);
    ll STRIDE=(argc>16)?atoll(argv[16]):1, SHIFT=(argc>17)?atoll(argv[17]):0;

    /* class consistency:  RD*(J1+J2) = (RD + 2*RPN*RD/RPD + RRN*RD/RRD)*M ; use rationals */
    /* 1 + 2 rho_p + rho_r = (J1+J2)/M   <=>  (RPD*RRD + 2*RPN*RRD + RRN*RPD)*M == RPD*RRD*(J1+J2) */
    if((RPD*RRD + 2*RPN*RRD + RRN*RPD)*M != RPD*RRD*(J1+J2)){
        fprintf(stderr,"# BAD CLASS: (J1+J2)/M != 1+2rho_p+rho_r\n"); return 2; }

    int prs[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,0};
    NP=0;
    for(int i=0;prs[i];i++){
        int p=prs[i]; if(p>N) break;
        int need = 2*vpfact(N,p)+2;
        u128 Mo=1; int K=0;
        while(K<need){ if(Mo > (((u128)1)<<126)/p) break; Mo*=(u128)p; K++; }
        PR[NP]=p; MOD[NP]=Mo;
        for(int n=0;n<=N+1;n++) VPF[NP][n]=vpfact(n,p);
        NP++;
    }
    NPS=0;
    for(int i=0;prs[i];i++){
        int p=prs[i]; if(p>NS) break;
        int need = 2*vpfact(NS,p)+2;
        unsigned long long Mo=1; int K=0;
        while(K<need && Mo <= (1ULL<<31)/p){ Mo*=(unsigned long long)p; K++; }
        if(K<need) continue;
        PRS[NPS]=p; MODS[NPS]=Mo;
        for(int n=0;n<=NS+1;n++) VPFS[NPS][n]=vpfact(n,p);
        NPS++;
    }

    ll M2 = M*M;
    ll hits=0, cnt=0;
    fprintf(stderr,"# mixed class rho_p=%lld/%lld rho_r=%lld/%lld M=%lld J=(%lld,%lld) r in [%lld,%lld] stride %lld/%lld\n",
            RPN,RPD,RRN,RRD,M,J1,J2,RMIN,RMAX,SHIFT,STRIDE);

    for(ll r=RMIN; r<=RMAX; r++){
        if(r==0) continue;
        if(((r-RMIN)%STRIDE) != SHIFT) continue;
        fflush(stdout); fprintf(stderr,"# r=%lld cnt=%lld hits=%lld\n",r,cnt,hits); fflush(stderr);
        /* p = C M^2 / r  must be an integer -> C multiple of Cstep = r/gcd(r,M^2) */
        ll Cstep = r/gcdll(r,M2); if(Cstep<0) Cstep=-Cstep;
        for(ll C=-GMAX; C<=GMAX; C++){
            if(C==0) continue;
            if(((C%Cstep)+Cstep)%Cstep) continue;
            ll g = C*M2;
            if(g % r) continue;
            ll pp = g/r;                       /* p = lam2*lam3 */
            ll hh = -(1 + 2*RPN*RRD/RPD*0) ;   /* placeholder, computed below */
            /* h = -(1+2rho_p+rho_r) g = -((J1+J2)/M) g */
            if((g*(J1+J2)) % M) continue;
            hh = -(g*(J1+J2))/M;
            ll jj;
            { __int128 t = (__int128)C*J1*J2; jj = (ll)t; }
            for(ll a=-AMAX; a<=AMAX; a++){
                ll s = a - r;                  /* s = lam2+lam3 */
                ll d = s*r + pp;
                if(d>DMAX || d<-DMAX) continue;
                /* b = (1-rho_r) r + (1-rho_p) s   ;  e = -rho_p(2p + r s) - rho_r r s */
                __int128 bn = (__int128)(RRD-RRN)*r*RPD + (__int128)(RPD-RPN)*s*RRD;  /* * RRD*RPD */
                if(bn % ((__int128)RRD*RPD)) continue;
                ll b = (ll)(bn/((__int128)RRD*RPD));
                __int128 en = -(__int128)RPN*RRD*(2*pp + r*s) - (__int128)RRN*RPD*r*s;
                if(en % ((__int128)RRD*RPD)) continue;
                ll e = (ll)(en/((__int128)RRD*RPD));
                /* discriminant of chi */
                {   __int128 D = (__int128)18*a*d*g - (__int128)4*a*a*a*g
                               + (__int128)a*a*d*d - (__int128)4*d*d*d - (__int128)27*g*g;
                    if(D==0) continue; }
                ll R2 = 4*g + 2*hh + jj;       /* R(2) */
                for(ll c=-CMAX; c<=CMAX; c++){
                    /* U_2 = P(1)c - Q(1) = (a+b+c)c - (d+e) - f  ; need 4 | U_2 */
                    ll K2 = (a+b+c)*c - (d+e);
                    ll r1,m1,r2,m2,rr,mm;
                    if(!lincong(-1, -K2, 4, &r1, &m1)) continue;   /* f == K2 (mod 4) */
                    /* U_3 = P(2)U_2 - 4 Q(2) c + 4 R(2),  Q(2)=4d+2e+f, U_2 = K2 - f */
                    ll P2c = 4*a+2*b+c;
                    ll K3 = P2c*K2 - 4*(4*d+2*e)*c + 4*R2;
                    ll A3 = -(P2c + 4*c);
                    if(!lincong(A3, -K3, 36, &r2, &m2)) continue;
                    if(!crt2(r1,m1,r2,m2,&rr,&mm)) continue;
                    ll fst = -FMAX;
                    { ll t=((fst-rr)%mm+mm)%mm; fst-=t; if(fst<-FMAX) fst+=mm; }
                    for(ll f=fst; f<=FMAX; f+=mm){
                        ll U2 = K2 - f;
                        ll U3 = P2c*U2 - 4*(4*d+2*e+f)*c + 4*R2;
                        if(U3 % 36) continue;
                        /* U_4 = P(3)U_3 - 9 Q(3) U_2 + 36 R(3) c */
                        __int128 U4 = (__int128)(9*a+3*b+c)*U3 - (__int128)9*(9*d+3*e+f)*U2
                                    + (__int128)36*(9*g+3*hh+jj)*c;
                        if(U4 % 576) continue;
                        if(c==0 && U2==0 && U3==0) continue;
                        if(a==c && b==2*a && d==f && e==2*d && hh==2*g && jj==g) continue;
                        if(!nontrivial(a,b,c,d,e,f,g,hh,jj)) continue;
                        if(!test_small(a,b,c,d,e,f,g,hh,jj)) continue;
                        cnt++;
                        int ok=1;
                        for(int i=0;i<NP;i++) if(!testp(a,b,c,d,e,f,g,hh,jj,i,N)){ok=0;break;}
                        if(ok){ hits++; printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld\n",
                                 RPN,RPD,RRN,RRD,M,J1,J2,r,a,c,d,f,C); }
                    }
                }
            }
        }
    }
    fprintf(stderr,"# deep %lld hits %lld\n",cnt,hits);
    return 0;
}
