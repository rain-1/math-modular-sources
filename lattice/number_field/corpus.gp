default(parisizemax, 4000000000);
default(realprecision, 60);
NK = 150;
LCT = vector(NK+2); LCT[1]=1; for(n=1,NK+1, LCT[n+1]=lcm(LCT[n],n));
\\ ---- generic (n+1)^r u_{n+1} = P(n) u_n - Q(n) u_{n-1} ------------------
seqrow(r,Pf,Qf,N) = my(v=vector(N+2)); v[1]=1; v[2]=Pf(0); for(n=1,N, v[n+2]=(Pf(n)*v[n+1]-Qf(n)*v[n])/(n+1)^r); v;
comprow(r,Pf,Qf,N) = my(v=vector(N+2)); v[1]=0; v[2]=1; for(n=1,N, v[n+2]=(Pf(n)*v[n+1]-Qf(n)*v[n])/(n+1)^r); v;
badat(b,k,N) = my(bad=-1); for(n=1,N, if(denominator(LCT[n+1]^k*b[n+1])!=1, bad=n; break)); bad;
sharpk(b,N) = my(kk=-1); for(k=0,8, if(badat(b,k,N)<0, kk=k; break)); kk;
introw(v,N) = my(bad=-1); for(n=0,N, if(denominator(v[n+1])!=1, bad=n; break)); bad;
\\ ---- reporting ---------------------------------------------------------
rep(nm, r, Pf, Qf, cp) = my(b,u,kk,rts,ab,l1,l2,sc,N=NK); u=seqrow(r,Pf,Qf,N); b=comprow(r,Pf,Qf,N); kk=sharpk(b,N); rts=polroots(cp); ab=vecsort(vector(#rts,i,abs(rts[i])),,4); l1=ab[1]; l2=if(#Set(ab)>1, vecsort(Set(ab),,4)[2], ab[1]); sc=-log(l2)-kk; print("=== ",nm); print("  charpoly            : ",cp); print("  factorisation over Q: ",factor(cp)); print("  disc                : ",poldisc(cp)); print("  roots               : ",rts); print("  |roots|             : ",ab); print("  distinct |roots|    : ",vecsort(Set(ab),,4)); print("  lambda_1, lambda_2  : ",l1,", ",l2); print("  row integral? first non-integral n<=",N," : ",introw(u,N)); print("  u_0..u_5            : ",vector(6,i,u[i])); print("  sharp k (companion b0=0,b1=1, n<=",N,") = ",kk,"   first fail at k-1: n=",badat(b,kk-1,N)); print("  rho = 1/|lambda_2|  : ",1/l2); print("  score = -log|l2| - k: ",sc); print("  complex pair?       : ",if(poldisc(cp)<0,"YES (|l1|=|l2|, no fold possible)","no")); print("");
\\ ------------------------- rows 1-7 (second order) -----------------------
rep("1. Zagier B  P=9n^2+9n+3, Q=27n^2", 2, (n)->9*n^2+9*n+3, (n)->27*n^2, x^2-9*x+27);
rep("2. Zagier D (Apery zeta(2))  P=11n^2+11n+3, Q=-n^2", 2, (n)->11*n^2+11*n+3, (n)->-n^2, x^2-11*x-1);
rep("3. Zagier E  P=12n^2+12n+4, Q=32n^2", 2, (n)->12*n^2+12*n+4, (n)->32*n^2, x^2-12*x+32);
rep("4. Herfurtner #4  P=117n^2+78n+21, Q=441(3n-1)^2", 2, (n)->117*n^2+78*n+21, (n)->441*(3*n-1)^2, x^2-117*x+3969);
rep("5. Herfurtner #6  P=72n^2+36n+6, Q=108(4n-1)(4n-3)", 2, (n)->72*n^2+36*n+6, (n)->108*(4*n-1)*(4*n-3), x^2-72*x+1728);
rep("6. sqrt(AZ(7,3,81))  P=56n^2+28n+6, Q=324(2n-1)^2", 2, (n)->56*n^2+28*n+6, (n)->324*(2*n-1)^2, x^2-56*x+1296);
rep("7. Beukers/sqrt(Apery)  P=136n^2+68n+10, Q=4(2n-1)^2", 2, (n)->136*n^2+68*n+10, (n)->4*(2*n-1)^2, x^2-136*x+16);
\\ ------------------------- rows 8-11 (cube type) -------------------------
rep("8. Apery zeta(3)", 3, (n)->(2*n+1)*(17*n^2+17*n+5), (n)->n^3, x^2-34*x+1);
rep("9. AZ eta (11,5,125)", 3, (n)->(2*n+1)*(11*n^2+11*n+5), (n)->125*n^3, x^2-22*x+125);
rep("10. AZ delta (7,3,81)", 3, (n)->(2*n+1)*(7*n^2+7*n+3), (n)->81*n^3, x^2-14*x+81);
rep("11. AZ(9,3,-27)", 3, (n)->(2*n+1)*(9*n^2+9*n+3), (n)->-27*n^3, x^2-18*x-27);
