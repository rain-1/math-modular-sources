/* rq_scan.c -- exhaustive scan for Apery-like rows over real quadratic fields.
   Shapes:
     R3: (n+1)^3 u_{n+1} = (2n+1)(A n^2+A n+B) u_n - C n^3 u_{n-1}
     R2: (n+1)^2 u_{n+1} = (A n^2+A n+B) u_n - C n^2 u_{n-1}
   Integrality test via v_n = (n!)^k u_n which satisfies the DIVISION-FREE recurrence
     R3: v_{n+1} = (2n+1)(A n^2+A n+B) v_n - C n^6 v_{n-1}
     R2: v_{n+1} =        (A n^2+A n+B) v_n - C n^4 v_{n-1}
   and u_n in O_K  <=>  (n!)^k | v_n in O_K.  All arithmetic modular & exact.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

typedef unsigned long long u64;
typedef unsigned __int128 u128;
typedef long long i64;

static int KEXP;       /* 3 for R3, 2 for R2 */
static int QPOW;       /* 6 for R3, 4 for R2  (exponent of n in Q(n)=C n^QPOW) */
static int TRFAC;      /* 1 -> multiply P by (2n+1); 0 -> no */

#define NMAXCHK 34

/* ---------- generic modular O_K arithmetic, modulus m (0 means 2^64) ------- */
static inline u64 mmul(u64 a,u64 b,u64 m){ if(!m) return a*b; return (u64)(((u128)a*b)%m); }
static inline u64 madd(u64 a,u64 b,u64 m){ if(!m) return a+b; u64 r=a+b; if(r>=m) r-=m; return r; }
static inline u64 msub(u64 a,u64 b,u64 m){ if(!m) return a-b; return a>=b? a-b : a+m-b; }

/* p-adic valuations of n! */
static int vpfact[64][32]; /* [pidx][n] */
static int primes[]={2,3,5,7,11,13,17,19,23,29};
#define NPR 10
static u64 pmod[NPR];   /* p^e modulus */
static int pexp[NPR];   /* e */
static u64 ppow[NPR][64];

static void initpv(void){
  for(int i=0;i<NPR;i++){
    int p=primes[i];
    /* largest e with p^e < 2^63 */
    u64 v=1; int e=0;
    while(v <= (u64)0x7fffffffffffffffULL / (u64)p){ v*=p; e++; }
    pmod[i]= (p==2)?0ULL:v;  /* p=2 handled by native wraparound mod 2^64 */
    pexp[i]= (p==2)?64:e;
    ppow[i][0]=1; for(int j=1;j<64;j++){ if(j<=pexp[i] && pexp[i]<64) ppow[i][j]=ppow[i][j-1]*p; else ppow[i][j]=0; }
    if(p==2){ for(int j=0;j<64;j++) ppow[i][j]= (j<64)?(1ULL<<j):0; }
    for(int n=0;n<32;n++){
      int tot=0; for(int q=2;q<=n;q++){ int x=q; while(x%p==0){tot++;x/=p;} }
      vpfact[i][n]=tot;
    }
  }
}

/* run the v-recurrence mod m and check p^{KEXP*vp(n!)} | v_n for n<=nmax.
   returns 1 if all pass. */
static inline int padic_ok(int pi,int nmax,u64 t,u64 s,
                           u64 a0,u64 a1,u64 b0,u64 b1,u64 c0,u64 c1)
{
  u64 m = pmod[pi]; int p = primes[pi]; int e = pexp[pi];
  u64 vp0=0,vp1=0, v0=1,v1=0;
  for(int n=0;n<nmax;n++){
    u64 nn = (u64)((i64)n*(i64)(n+1));
    if(m) nn %= m;
    /* P = A*n(n+1) + B ; times (2n+1) if TRFAC */
    u64 P0 = madd(mmul(a0,nn,m), b0, m);
    u64 P1 = madd(mmul(a1,nn,m), b1, m);
    if(TRFAC){ u64 f=(u64)(2*n+1); if(m) f%=m; P0=mmul(P0,f,m); P1=mmul(P1,f,m); }
    /* Q = C * n^QPOW */
    u64 nq=1,nb=(u64)n; if(m) nb%=m;
    for(int j=0;j<QPOW;j++) nq=mmul(nq,nb,m);
    u64 Q0=mmul(c0,nq,m), Q1=mmul(c1,nq,m);
    /* new = P*v - Q*vp   (O_K mult: (x0,x1)(y0,y1)=(x0y0+s x1y1, x0y1+x1y0+t x1y1)) */
    u64 r0,r1,g0,g1;
    r0 = madd(mmul(P0,v0,m), mmul(s,mmul(P1,v1,m),m), m);
    r1 = madd(madd(mmul(P0,v1,m),mmul(P1,v0,m),m), mmul(t,mmul(P1,v1,m),m), m);
    g0 = madd(mmul(Q0,vp0,m), mmul(s,mmul(Q1,vp1,m),m), m);
    g1 = madd(madd(mmul(Q0,vp1,m),mmul(Q1,vp0,m),m), mmul(t,mmul(Q1,vp1,m),m), m);
    vp0=v0; vp1=v1;
    v0=msub(r0,g0,m); v1=msub(r1,g1,m);
    int need = KEXP*vpfact[pi][n+1];
    if(need>0){
      if(need>e) return 1;            /* cannot certify further: stop, accept */
      if(p==2){ u64 msk=(need>=64)?~0ULL:((1ULL<<need)-1); if((v0&msk)||(v1&msk)) return 0; }
      else { u64 q=ppow[pi][need]; if(v0%q || v1%q) return 0; }
    }
  }
  return 1;
}

/* ------------------------- sieve tables --------------------------------- */
typedef struct { int m; int nel; int stride; unsigned short *lst; int *cnt; } Sieve;

static void build_sieve(Sieve*S,int m,int pi,int nmax,u64 t,u64 s,u64 c0,u64 c1){
  int nel=m*m; S->m=m; S->nel=nel; S->stride=nel;
  S->lst=(unsigned short*)malloc((size_t)nel*nel*sizeof(unsigned short));
  S->cnt=(int*)calloc(nel,sizeof(int));
  u64 mm=(u64)m; u64 tm=t%mm, sm=s%mm, cc0=c0%mm, cc1=c1%mm;
  int p=primes[pi];
  #pragma omp parallel for schedule(dynamic,8)
  for(int ia=0;ia<nel;ia++){
    u64 a0=ia/m, a1=ia%m; int k=0;
    for(int ib=0;ib<nel;ib++){
      u64 b0=ib/m, b1=ib%m;
      /* inline small run */
      u64 vp0=0,vp1=0,v0=1,v1=0; int ok=1;
      for(int n=0;n<nmax;n++){
        u64 nn=((u64)((i64)n*(i64)(n+1)))%mm;
        u64 P0=(a0*nn+b0)%mm, P1=(a1*nn+b1)%mm;
        if(TRFAC){ u64 f=(u64)(2*n+1)%mm; P0=P0*f%mm; P1=P1*f%mm; }
        u64 nq=1,nb=((u64)n)%mm; for(int j=0;j<QPOW;j++) nq=nq*nb%mm;
        u64 Q0=cc0*nq%mm, Q1=cc1*nq%mm;
        u64 r0=(P0*v0 + sm*((P1*v1)%mm))%mm;
        u64 r1=(P0*v1 + P1*v0 + tm*((P1*v1)%mm))%mm;
        u64 g0=(Q0*vp0 + sm*((Q1*vp1)%mm))%mm;
        u64 g1=(Q0*vp1 + Q1*vp0 + tm*((Q1*vp1)%mm))%mm;
        vp0=v0; vp1=v1;
        v0=(r0+mm-g0)%mm; v1=(r1+mm-g1)%mm;
        int need=KEXP*vpfact[pi][n+1];
        if(need>0){ u64 q=1; for(int j=0;j<need;j++) q*=p;
                    if(q<=mm){ if(v0%q||v1%q){ ok=0; break; } } }
      }
      if(ok) S->lst[(size_t)ia*nel + (k++)]=(unsigned short)ib;
    }
    S->cnt[ia]=k;
  }
}
static void free_sieve(Sieve*S){ free(S->lst); free(S->cnt); }

/* ---------------------- decimal string mod m ----------------------------- */
static u64 strmod(const char*str,u64 m){
  int neg=0; const char*q=str; if(*q=='-'){neg=1;q++;} if(*q=='+')q++;
  u64 r=0;
  for(;*q;q++){ if(*q<'0'||*q>'9') break;
    if(m) r=(u64)(((u128)r*10 + (u64)(*q-'0'))%m); else r=r*10+(u64)(*q-'0'); }
  if(neg) r = m? (m-r)%m : (u64)(0-r);
  return r;
}

/* ------------------------------ main ------------------------------------- */
typedef struct { int D,t,s; } Fld;

int main(int argc,char**argv){
  if(argc<7){ fprintf(stderr,"usage: rq_scan shape(2|3) fieldsfile Cmark(C|S) ABOX BBOX out\n"); return 1; }
  int shape=atoi(argv[1]);
  const char*ff=argv[2]; char cmark=argv[3][0];
  int ABOX=atoi(argv[4]), BBOX=atoi(argv[5]);
  FILE*out=fopen(argv[6],"w");
  if(shape==3){ KEXP=3; QPOW=6; TRFAC=1; } else { KEXP=2; QPOW=4; TRFAC=0; }
  initpv();
  int M2 = (shape==3)?8:4;
  int M3 = (shape==3)?27:9;
  int NM2 = (shape==3)?3:3;   /* nmax certifiable mod M2 */
  int NM3 = (shape==3)?5:5;   /* nmax certifiable mod M3 */
  int nmax2 = (shape==3)?23:33;
  int nmax3 = (shape==3)?29:30;
  int nmaxo = 30;

  /* read fields file */
  FILE*f=fopen(ff,"r"); if(!f){perror("fields");return 1;}
  char line[4096];
  int curD=0,curt=0,curs=0;
  char (*Cs)[2][256]=NULL; long double *Cv1=NULL,*Cv2=NULL; int nC=0, capC=0;
  long double MITHR = (argc>7)? strtold(argv[7],NULL) : -1e30L;
  int KFIL = (argc>8)? atoi(argv[8]) : ((shape==3)?3:2);
  long long tested=0, integral=0;
  while(fgets(line,sizeof line,f)){
    if(!strncmp(line,"D ",2)){
      curD=0;curt=0;curs=0; sscanf(line,"D %d %d %d",&curD,&curt,&curs);
      nC=0;
    } else if(!strncmp(line,"NC1 ",4)||!strncmp(line,"NC2 ",4)){
      /* nothing */
    } else if(line[0]==cmark && line[1]==' '){
      if(nC>=capC){ capC=capC?capC*2:64; Cs=realloc(Cs,(size_t)capC*sizeof(*Cs));
        Cv1=realloc(Cv1,(size_t)capC*sizeof(long double)); Cv2=realloc(Cv2,(size_t)capC*sizeof(long double)); }
      char x[256],y[256],e1[256],e2[256];
      int nf=sscanf(line,"%*s %255s %255s %255s %255s",x,y,e1,e2);
      strcpy(Cs[nC][0],x); strcpy(Cs[nC][1],y);
      Cv1[nC]= (nf>=3)? strtold(e1,NULL):0.0L; Cv2[nC]= (nf>=4)? strtold(e2,NULL):0.0L; nC++;
    } else if(!strncmp(line,"ENDD",4) || (!strncmp(line,"NC2 ",4)&&cmark=='C')){
      /* handled below by flush at next D — see flush logic */
    }
    /* flush when we hit the marker that ends this D's relevant C list */
    int doflush = (!strncmp(line,"ENDD",4) && nC>0);
    if(!doflush) continue;

    /* ---- scan this D ---- */
    int D=curD; u64 t=(u64)curt, s=(u64)curs;
    long double WV1=((long double)curt+sqrtl((long double)(curt*curt+4*curs)))/2.0L;
    long double WV2=((long double)curt-sqrtl((long double)(curt*curt+4*curs)))/2.0L;
    for(int ci=0;ci<nC;ci++){
      u64 c2_0=strmod(Cs[ci][0],0), c2_1=strmod(Cs[ci][1],0);   /* mod 2^64 */
      u64 cm0[NPR],cm1[NPR];
      for(int i=0;i<NPR;i++){ cm0[i]=strmod(Cs[ci][0],pmod[i]); cm1[i]=strmod(Cs[ci][1],pmod[i]); }
      long double CV1=Cv1[ci], CV2=Cv2[ci];
      Sieve S2,S3;
      build_sieve(&S2,M2,0,NM2,t,s,strmod(Cs[ci][0],(u64)M2),strmod(Cs[ci][1],(u64)M2));
      build_sieve(&S3,M3,1,NM3,t,s,strmod(Cs[ci][0],(u64)M3),strmod(Cs[ci][1],(u64)M3));
      int MM = M2*M3;
      /* CRT helper table: for r2 in [0,M2), r3 in [0,M3) -> r mod MM */
      static int crt[64][64];
      for(int r2=0;r2<M2;r2++) for(int r3=0;r3<M3;r3++){
        int x; for(x=0;x<MM;x++) if(x%M2==r2 && x%M3==r3) break; crt[r2][r3]=x; }
      long long ltested=0, lintegral=0;
      #pragma omp parallel for schedule(dynamic,8) reduction(+:ltested,lintegral)
      for(int a0=-ABOX;a0<=ABOX;a0++){
        char buf[8192]; int blen=0;
        for(int a1=-ABOX;a1<=ABOX;a1++){
          if(MITHR > -1e29L){
            long double bestm=-1e30L;
            for(int pl=0;pl<2;pl++){
              long double aa[2],cc[2],L1[2],L2[2];
              aa[0]=(long double)a0+(long double)a1*WV1; aa[1]=(long double)a0+(long double)a1*WV2;
              cc[0]=CV1; cc[1]=CV2;
              for(int q=0;q<2;q++){
                long double A_=(shape==3)?aa[q]:aa[q]/2.0L, c_=cc[q], d_=A_*A_-c_;
                if(d_>0){ long double r1=fabsl(A_)+sqrtl(d_); L1[q]=logl(r1); L2[q]=logl(fabsl(c_))-L1[q]; }
                else { L1[q]=0.5L*logl(fabsl(c_)); L2[q]=L1[q]; }
              }
              int vo=1-pl;
              long double m = -(L2[pl]+L1[vo])/2.0L - (long double)KFIL;
              if(m>bestm) bestm=m;
            }
            if(bestm < MITHR) continue;
          }
          ltested += (long long)(2*BBOX+1)*(2*BBOX+1);
          int ia2 = ((a0%M2+M2)%M2)*M2 + ((a1%M2+M2)%M2);
          int ia3 = ((a0%M3+M3)%M3)*M3 + ((a1%M3+M3)%M3);
          int n2=S2.cnt[ia2], n3=S3.cnt[ia3];
          if(!n2||!n3) continue;
          for(int i2=0;i2<n2;i2++){
            int ib2=S2.lst[(size_t)ia2*S2.nel+i2]; int r20=ib2/M2, r21=ib2%M2;
            for(int i3=0;i3<n3;i3++){
              int ib3=S3.lst[(size_t)ia3*S3.nel+i3]; int r30=ib3/M3, r31=ib3%M3;
              int R0=crt[r20][r30], R1=crt[r21][r31];
              /* enumerate b0 = R0 + k*MM in [-BBOX,BBOX] */
              int st0 = R0 - MM*((R0+BBOX)/MM);
              for(int b0=st0;b0<=BBOX;b0+=MM){ if(b0< -BBOX) continue;
                int st1 = R1 - MM*((R1+BBOX)/MM);
                for(int b1=st1;b1<=BBOX;b1+=MM){ if(b1< -BBOX) continue;
                  u64 A0=(u64)(i64)a0, A1=(u64)(i64)a1, B0=(u64)(i64)b0, B1=(u64)(i64)b1;
                  if(!padic_ok(0,nmax2,t,s,A0,A1,B0,B1,c2_0,c2_1)) continue;
                  int good=1;
                  for(int pi=1;pi<NPR && good;pi++){
                    u64 mm=pmod[pi];
                    u64 aa0=(u64)(((i64)a0%(i64)mm+(i64)mm)%(i64)mm);
                    u64 aa1=(u64)(((i64)a1%(i64)mm+(i64)mm)%(i64)mm);
                    u64 bb0=(u64)(((i64)b0%(i64)mm+(i64)mm)%(i64)mm);
                    u64 bb1=(u64)(((i64)b1%(i64)mm+(i64)mm)%(i64)mm);
                    int nx = (pi==1)?nmax3:nmaxo;
                    if(!padic_ok(pi,nx,t%mm,s%mm,aa0,aa1,bb0,bb1,cm0[pi],cm1[pi])) good=0;
                  }
                  if(good){ lintegral++;
                    blen+=snprintf(buf+blen,sizeof(buf)-blen,"R%d %d %d %d %d %d %s %s\n",
                                   shape,D,a0,a1,b0,b1,Cs[ci][0],Cs[ci][1]);
                    if(blen>7000){ 
                      #pragma omp critical
                      { fputs(buf,out); fflush(out);} blen=0; buf[0]=0; }
                  }
                }
              }
            }
          }
        }
        if(blen){ 
          #pragma omp critical
          { fputs(buf,out); fflush(out);} }
      }
      tested+=ltested; integral+=lintegral;
      free_sieve(&S2); free_sieve(&S3);
    }
    fprintf(stderr,"D=%d done  tested=%lld integral=%lld\n",D,tested,integral);
    nC=0;
  }
  fprintf(stderr,"TOTAL tested=%lld integral=%lld\n",tested,integral);
  printf("TOTAL tested=%lld integral=%lld\n",tested,integral);
  fclose(out);
  return 0;
}
