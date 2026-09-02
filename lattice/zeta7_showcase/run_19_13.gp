default(parisize,"24G");
read("/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/lift_lib.gp");
RR=19;
DD=13;
TAG="19_13";
NPMAX=300;
TP=(L+1)^7*(L^2-14*L+1)^6;
tv=vector(20,i,polcoeff(TP,20-i));
EXROWS=matrix(19,280,a,b,my(v=0);if(b==(a+1)*14,v=1);if(b==14,v=v-tv[a+1]);v);
ps=mkprimes(NPMAX);
ACC=0;
piv0=0;
QH=vector(NPMAX);
OKH=vector(NPMAX);
nstab=0;
gettime();
{for(i=1,NPMAX, setC(ps[i]); res=canon(RR,DD,ps[i]); if(i==1, piv0=res[2]; ACC=res[1]*Mod(1,ps[i]), if(res[2]!=piv0, error("pivot mismatch at prime ",i)); ACC=chinese(ACC,res[1]*Mod(1,ps[i]))); Q=bestappr(ACC); OKH[i]=(type(Q)=="t_MAT"); if(OKH[i],QH[i]=Q); if(i>=4 && OKH[i] && OKH[i-3], if(QH[i]==QH[i-3], nstab=i; break)));}
print("pivots: ",piv0);
print("stable at prime count = ",nstab,"   time ",gettime()," ms");
if(nstab==0,print("NOT STABLE with ",NPMAX," primes"));
QF=QH[nstab];
KDIM=matsize(QF)[1];
print("kernel dim = ",KDIM);
BV=vector(KDIM,i,prim(QF[i,]));
{for(i=1,KDIM, print("basis ",i," max coeff digits = ",maxdig(BV[i])));}
PVL=vector(KDIM,i,tovec2pol(BV[i],RR,DD));
{for(i=1,KDIM, print("basis ",i," EXACT verify residual count (should be 0) = ",verifyrec(PVL[i],RR)));}
write(concat(BASE,"rec_19_13_raw.gp"),Str("BVL=",BV));
print("done");
quit;
