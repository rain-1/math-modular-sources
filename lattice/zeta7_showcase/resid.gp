default(parisize,"24G");
default(linewrap,0);
BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
DN=read(concat(BASE,"dn.txt"));
ND=400;
resid(PV,r,nmax)=vector(nmax-r,i,my(nn=r+i,s=0);for(j=0,r,s=s+subst(PV[j+1],n,nn)*DN[nn-j+1]);s);
read(concat(BASE,"rec_19_13.gp"));
gettime();
R=resid(PV,19,400);
print("(19,13): computed ",#R," residuals, t=",gettime());
print("first 6 R_n (n=20..25):");
{for(i=1,6,print("  n=",19+i,"  num-digits=",#Str(abs(numerator(R[i])))," den-digits=",#Str(abs(denominator(R[i])))," sign=",sign(R[i])));}
print("R_21/R_20 = ",R[2]/R[1]);
print("R_22/R_20 = ",R[3]/R[1]);
print("all |R_n| equal? ",vecmax(vector(#R,i,abs(R[i])!=abs(R[1])))==0);
print("R_n = R_20*(-1)^n pattern? ",vecmax(vector(#R,i,R[i]!=(-1)^(i-1)*R[1]))==0);
print("R_n constant? ",vecmax(vector(#R,i,R[i]!=R[1]))==0);
print("R_20 = ",R[1]);
quit;
