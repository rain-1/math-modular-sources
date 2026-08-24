/* FIVE-term (rank 2, SIX singular point) integrality scan.
 *
 * Row: (n+1)^2 u_{n+1} = P(n)u_n - Q(n)u_{n-1} + R(n)u_{n-2} - T(n)u_{n-3},
 *      u_0=1, u_{-1}=u_{-2}=u_{-3}=0,
 *      P=a n^2+b n+c, Q=d n^2+e n+f, R=g n^2+h n+j, T=k n^2+l n+m.
 * L = theta^2 - tP(theta) + t^2 Q(theta+1) - t^3 R(theta+2) + t^4 T(theta+3)
 *   = t^2 Rc D^2 + t Sc D + t Vc,
 *   Rc = 1 - a t + d t^2 - g t^3 + k t^4
 *   Sc = 1 - (a+b)t + (3d+e)t^2 - (5g+h)t^3 + (7k+l)t^4
 *   Vc = -c + (d+e+f)t - (4g+2h+j)t^2 + (9k+3l+m)t^3
 * Exponents (0,0) at t=0, (0,rho_i) at the four roots of Rc, (3-s_1,3-s_2) at
 * infinity with s_i the roots of T.  Equal rho at all four finite points iff
 *   b=(1-rho)a, e=-2 rho d, h=-(1+3rho)g, l=-(2+4rho)k,
 * and Fuchs reads sum rho_i = s_1+s_2-2, i.e. (J1+J2)/M = 2+4 rho for
 *   T(n) = C(Mn-J1)(Mn-J2),  k=C M^2, l=-C M(J1+J2), m=C J1 J2.
 * Free parameters: a,c,d,f,g,j,C  (three accessory parameters c,f,j -- as a
 * rank-2 Fuchsian equation with six singular points must have).
 *
 * Integrality tested on U_n = u_n (n!)^2:
 *   U_{n+1} = P U_n - n^2 Q U_{n-1} + n^2(n-1)^2 R U_{n-2}
 *                                   - n^2(n-1)^2(n-2)^2 T U_{n-3}.
 *
 * usage: 04_fscan5 RN RD M J1 J2 AMIN AMAX CMAX DMAX FMAX GMAX JMAX KMAX N [STRIDE SHIFT]
 * out:   RN RD M J1 J2 a c d f g j C
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

#define MP1 2305843009213693951ULL
#define MP2 2305843009213693739ULL
static inline unsigned long long mmul(unsigned long long x,unsigned long long y,unsigned long long m)
{ return (unsigned long long)(((u128)x*y)%m); }
static inline unsigned long long mred(ll x, unsigned long long m)
{ ll r = x % (ll)m; return (unsigned long long)(r<0? r+(ll)m : r); }
static unsigned long long mpow(unsigned long long b,unsigned long long e,unsigned long long m)
{ unsigned long long r=1; while(e){ if(e&1) r=mmul(r,b,m); b=mmul(b,b,m); e>>=1;} return r; }

/* discriminant of chi = lam^4 - a lam^3 + d lam^2 - g lam + k  (monic quartic
 * with (b,c,d,e) = (-a, d, -g, k) in the standard formula).  Zero iff chi has a
 * repeated root.  All terms fit comfortably in __int128 for the boxes used. */
static inline int quartic_disc_zero(ll A,ll D,ll G,ll K)
{
    __int128 b=-A, c=D, d=-G, e=K;
    __int128 b2=b*b, c2=c*c, d2=d*d, e2=e*e;
    __int128 T = (__int128)256*e*e2 - (__int128)192*b*d*e2 - (__int128)128*c2*e2
      + (__int128)144*c*d2*e - (__int128)27*d2*d2 + (__int128)144*b2*c*e2
      - (__int128)6*b2*d2*e - (__int128)80*b*c2*d*e + (__int128)18*b*c*d*d2
      + (__int128)16*c2*c2*e - (__int128)4*c2*c*d2 - (__int128)27*b2*b2*e2
      + (__int128)18*b2*b*c*d*e - (__int128)4*b2*b*d*d2 - (__int128)4*b2*c2*c*e
      + b2*c2*d2;
    return T==0;
}

static inline int test_small(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io,ll ka,ll la,ll mu)
{
    for(int i=0;i<NPS;i++){
        unsigned long long MM = MODS[i]; int p = PRS[i];
        unsigned long long r0=0,r1=0,r2=1,r3;   /* U_{-3},U_{-2},U_{-1}? */
        /* state: v[n-3],v[n-2],v[n-1],v[n] ; start n=0: U_{-3}=U_{-2}=U_{-1}=0,U_0=1 */
        unsigned long long w[4]={0,0,0,1};
        for(int n=0;n<=NS-1;n++){
            ll nn=n;
            ll Pn=(al*nn*nn+be*nn+ga)%(ll)MM, Qn=(de*nn*nn+ep*nn+ze)%(ll)MM;
            ll Rn=(et*nn*nn+th*nn+io)%(ll)MM, Tn=(ka*nn*nn+la*nn+mu)%(ll)MM;
            unsigned long long pP=(unsigned long long)(Pn<0?Pn+(ll)MM:Pn);
            unsigned long long pQ=(unsigned long long)(Qn<0?Qn+(ll)MM:Qn);
            unsigned long long pR=(unsigned long long)(Rn<0?Rn+(ll)MM:Rn);
            unsigned long long pT=(unsigned long long)(Tn<0?Tn+(ll)MM:Tn);
            unsigned long long c1=(unsigned long long)(nn*nn)%MM;
            unsigned long long c2=(unsigned long long)(nn*(nn-1))%MM; c2=(c2*c2)%MM;
            unsigned long long c3=(unsigned long long)(nn*(nn-1))%MM; c3=(c3*c3)%MM;
            c3=(c3*(((unsigned long long)((nn-2)*(nn-2)))%MM))%MM;
            unsigned long long nx=(pP*w[3])%MM;
            nx=(nx+MM-(((c1*pQ)%MM)*w[2])%MM)%MM;
            nx=(nx+(((c2*pR)%MM)*w[1])%MM)%MM;
            nx=(nx+MM-(((c3*pT)%MM)*w[0])%MM)%MM;
            int need=2*VPFS[i][n+1];
            if(need>0 && nx){
                unsigned long long x=nx; int v=0;
                while(v<need && x%(unsigned long long)p==0){x/=(unsigned long long)p;v++;}
                if(v<need) return 0;
            }
            w[0]=w[1];w[1]=w[2];w[2]=w[3];w[3]=nx;
        }
        (void)r0;(void)r1;(void)r2;(void)r3;
    }
    return 1;
}

static inline int testp(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io,ll ka,ll la,ll mu,int pi,int N)
{
    u128 MM=MOD[pi]; int p=PR[pi];
    #define RED(x) ( (x)>=0 ? (u128)(x)%MM : MM - ((u128)(-(x))%MM) )
    u128 w[4]={0,0,0,1};
    for(int n=0;n<=N-1;n++){
        ll nn=n;
        u128 pP=RED(al*nn*nn+be*nn+ga), pQ=RED(de*nn*nn+ep*nn+ze);
        u128 pR=RED(et*nn*nn+th*nn+io), pT=RED(ka*nn*nn+la*nn+mu);
        u128 c1=(u128)(nn*nn)%MM;
        u128 c2=(u128)(nn*(nn-1))%MM; c2=(c2*c2)%MM;
        u128 c3=(c2*((u128)((nn-2)*(nn-2))%MM))%MM;
        u128 nx=(pP*w[3])%MM;
        nx=(nx+MM-((((c1*pQ)%MM)*w[2])%MM))%MM;
        nx=(nx+((((c2*pR)%MM)*w[1])%MM))%MM;
        nx=(nx+MM-((((c3*pT)%MM)*w[0])%MM))%MM;
        int need=2*VPF[pi][n+1];
        if(need>0 && nx){
            u128 x=nx; int v=0;
            while(v<need && x%(u128)p==0){x/=(u128)p;v++;}
            if(v<need) return 0;
        }
        w[0]=w[1];w[1]=w[2];w[2]=w[3];w[3]=nx;
    }
    return 1;
    #undef RED
}

static inline int nontrivial5(ll al,ll be,ll ga,ll de,ll ep,ll ze,ll et,ll th,ll io,ll ka,ll la,ll mu)
{
    unsigned long long MM[2]={MP1,MP2}; int zr[2]={0,0};
    for(int q=0;q<2;q++){
        unsigned long long m=MM[q], w[4]={0,0,0,1}; int z=0;
        for(int n=0;n<=9;n++){
            unsigned long long pP=mred(al*(ll)n*n+be*n+ga,m), pQ=mred(de*(ll)n*n+ep*n+ze,m);
            unsigned long long pR=mred(et*(ll)n*n+th*n+io,m), pT=mred(ka*(ll)n*n+la*n+mu,m);
            unsigned long long c1=mred((ll)n*n,m);
            unsigned long long c2=mred((ll)n*(n-1),m); c2=mmul(c2,c2,m);
            unsigned long long c3=mmul(c2, mred((ll)(n-2)*(n-2),m), m);
            unsigned long long nx=mmul(pP,w[3],m);
            nx=(nx+m-mmul(mmul(c1,pQ,m),w[2],m))%m;
            nx=(nx+mmul(mmul(c2,pR,m),w[1],m))%m;
            nx=(nx+m-mmul(mmul(c3,pT,m),w[0],m))%m;
            if(nx==0){ if(++z>=4){zr[q]=1;break;} } else z=0;
            w[0]=w[1];w[1]=w[2];w[2]=w[3];w[3]=nx;
        }
    }
    return !(zr[0]&&zr[1]);
}

static ll modinv(ll a, ll m){
    if(m==1) return 0;
    ll g=m, xx=0, x1=1, a1=((a%m)+m)%m;
    while(a1){ ll q=g/a1,t=g-q*a1; g=a1;a1=t; t=xx-q*x1; xx=x1; x1=t; }
    return ((xx%m)+m)%m;
}
static inline int lincong(ll A, ll B, ll q, ll *r, ll *m){
    A=((A%q)+q)%q; B=((B%q)+q)%q;
    ll g=gcdll(A,q); if(g==0){ if(B) return 0; *r=0;*m=1;return 1; }
    if(B%g) return 0;
    ll q1=q/g; if(q1==1){*r=0;*m=1;return 1;}
    ll inv=modinv((A/g)%q1,q1);
    *r=(ll)(((__int128)((B/g)%q1)*inv)%q1); *m=q1; return 1;
}
static inline int crt2(ll r1,ll m1,ll r2,ll m2,ll *r,ll *m){
    ll g=gcdll(m1,m2); if(((r2-r1)%g+g)%g) return 0;
    ll m2g=m2/g, l=m1*m2g;
    ll inv=modinv((m1/g)%m2g,m2g);
    ll k=(ll)(((__int128)((((r2-r1)/g)%m2g+m2g)%m2g)*inv)%m2g);
    *r=(ll)((((__int128)r1+(__int128)m1*k)%l+l)%l); *m=l; return 1;
}

int main(int argc,char**argv)
{
    if(argc<15){fprintf(stderr,"usage: %s RN RD M J1 J2 AMIN AMAX CMAX DMAX FMAX GMAX JMAX KMAX N [STRIDE SHIFT]\n",argv[0]);return 1;}
    ll RN=atoll(argv[1]),RD=atoll(argv[2]),M=atoll(argv[3]),J1=atoll(argv[4]),J2=atoll(argv[5]);
    ll AMIN=atoll(argv[6]),AMAX=atoll(argv[7]),CMAX=atoll(argv[8]),DMAX=atoll(argv[9]);
    ll FMAX=atoll(argv[10]),GMAX=atoll(argv[11]),JMAX=atoll(argv[12]),KMAX=atoll(argv[13]);
    int N=atoi(argv[14]);
    ll STRIDE=(argc>15)?atoll(argv[15]):1, SHIFT=(argc>16)?atoll(argv[16]):0;
    int NODISC=(argc>17)?atoi(argv[17]):0;   /* 1: keep repeated-root rows too */

    if(RD*(J1+J2) != (2*RD+4*RN)*M){fprintf(stderr,"# BAD CLASS: (J1+J2)/M != 2+4rho\n");return 2;}

    int prs[]={2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113,0};
    NP=0;
    for(int i=0;prs[i];i++){int p=prs[i]; if(p>N)break;
        int need=2*vpfact(N,p)+2; u128 Mo=1; int K=0;
        while(K<need){ if(Mo > (((u128)1)<<126)/p) break; Mo*=(u128)p; K++; }
        PR[NP]=p;MOD[NP]=Mo; for(int n=0;n<=N+1;n++) VPF[NP][n]=vpfact(n,p); NP++;}
    NPS=0;
    for(int i=0;prs[i];i++){int p=prs[i]; if(p>NS)break;
        int need=2*vpfact(NS,p)+2; unsigned long long Mo=1; int K=0;
        while(K<need && Mo<= (1ULL<<31)/p){Mo*=(unsigned long long)p;K++;}
        if(K<need) continue;
        PRS[NPS]=p;MODS[NPS]=Mo; for(int n=0;n<=NS+1;n++) VPFS[NPS][n]=vpfact(n,p); NPS++;}

    ll Astep = RD/gcdll(RD, RD-RN);
    ll Dstep = RD/gcdll(RD, 2*RN);
    ll Gstep = RD/gcdll(RD, RD+3*RN);
    ll a0=AMIN; if(a0%Astep) a0 += Astep - ((a0%Astep)+Astep)%Astep;
    fprintf(stderr,"# 5-term class rho=%lld/%lld M=%lld J=(%lld,%lld) Astep=%lld Dstep=%lld Gstep=%lld\n",
            RN,RD,M,J1,J2,Astep,Dstep,Gstep);

    ll hits=0,cnt=0;
    ll W3=(3*M-J1)*(3*M-J2), W4=(4*M-J1)*(4*M-J2);
    for(ll a=a0;a<=AMAX;a+=Astep){
        if(((a-a0)/Astep)%STRIDE!=SHIFT) continue;
        fflush(stdout); fprintf(stderr,"# a=%lld cnt=%lld hits=%lld\n",a,cnt,hits); fflush(stderr);
        ll be=a*(RD-RN)/RD;
        for(ll c=-CMAX;c<=CMAX;c++){
            ll P1=a+be+c,P2=4*a+2*be+c,P3=9*a+3*be+c,P4=16*a+4*be+c;
            ll d0=-DMAX; if(d0%Dstep) d0 += Dstep - ((d0%Dstep)+Dstep)%Dstep;
            for(ll d=d0;d<=DMAX;d+=Dstep){
                ll ep=-2*RN*d/RD;
                for(ll f=-FMAX;f<=FMAX;f++){
                    ll Q1=d+ep+f,Q2=4*d+2*ep+f,Q3=9*d+3*ep+f,Q4=16*d+4*ep+f;
                    ll U2=P1*c-Q1; if(U2&3) continue;
                    ll g0=-GMAX; if(g0%Gstep) g0 += Gstep - ((g0%Gstep)+Gstep)%Gstep;
                    for(ll g=g0;g<=GMAX;g+=Gstep){
                        if(g==0) continue;
                        ll th=-(RD+3*RN)*g/RD;
                        /* U_3 = P(2)U_2 - 4Q(2)c + 4 R(2),  R(2)=4g+2th+j : linear in j */
                        ll K3=P2*U2-4*Q2*c+4*(4*g+2*th);
                        ll rj,mj;
                        if(!lincong(4,-K3,36,&rj,&mj)) continue;
                        ll jst=-JMAX; { ll t=((jst-rj)%mj+mj)%mj; jst-=t; if(jst<-JMAX) jst+=mj; }
                        for(ll j=jst;j<=JMAX;j+=mj){
                            ll R2=4*g+2*th+j, R3=9*g+3*th+j, R4=16*g+4*th+j;
                            ll U3=P2*U2-4*Q2*c+4*R2;
                            if(U3%36) continue;
                            /* U_4 = P(3)U_3 - 9Q(3)U_2 + 36 R(3) c - 36 T(3)   (T(3)=C W3) */
                            ll K4=P3*U3-9*Q3*U2+36*R3*c;  ll A4=-36*W3;
                            ll r1,m1,r2,m2,rr,mm;
                            if(!lincong(A4,-K4,576,&r1,&m1)) continue;
                            /* U_5 = P(4)U_4 - 16Q(4)U_3 + 144 R(4) U_2 - 576 T(4) U_1,
                               U_1=c, T(4)=C W4 ; linear in C.  need 2^6 3^2 5^2 = 14400 | U_5 */
                            ll K5=P4*K4-16*Q4*U3+144*R4*U2-0; ll A5=P4*A4-576*W4*c;
                            if(!lincong(A5,-K5,14400,&r2,&m2)) continue;
                            if(!crt2(r1,m1,r2,m2,&rr,&mm)) continue;
                            ll Cst=-KMAX; { ll t=((Cst-rr)%mm+mm)%mm; Cst-=t; if(Cst<-KMAX) Cst+=mm; }
                            for(ll C=Cst;C<=KMAX;C+=mm){
                                if(C==0) continue;
                                ll ka=C*M*M, la=-C*M*(J1+J2), mu=C*J1*J2;
                                /* constant-coefficient row */
                                if(a==c&&be==2*a&&d==f&&ep==2*d&&g==j&&th==2*g&&ka==mu&&la==2*ka) continue;
                                /* repeated root of the quartic */
                                int degen = quartic_disc_zero(a,d,g,ka);
                                if(degen && !NODISC) continue;
                                if(!nontrivial5(a,be,c,d,ep,f,g,th,j,ka,la,mu)) continue;
                                if(!test_small(a,be,c,d,ep,f,g,th,j,ka,la,mu)) continue;
                                cnt++;
                                int ok=1;
                                for(int i=0;i<NP;i++) if(!testp(a,be,c,d,ep,f,g,th,j,ka,la,mu,i,N)){ok=0;break;}
                                if(ok){hits++;printf("%lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld %lld%s\n",
                                        RN,RD,M,J1,J2,a,c,d,f,g,j,C, degen?"  # REPEATED":"");}
                            }
                        }
                    }
                }
            }
        }
    }
    fprintf(stderr,"# deep %lld hits %lld\n",cnt,hits);
    return 0;
}
