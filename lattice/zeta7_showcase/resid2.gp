default(parisize,"24G");
default(linewrap,0);
BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
DN=read(concat(BASE,"dn.txt"));
resid(PV,r,nmax)=vector(nmax-r,i,my(nn=r+i,s=0);for(j=0,r,s=s+subst(PV[j+1],n,nn)*DN[nn-j+1]);s);
read(concat(BASE,"rec_19_13.gp"));
R=resid(PV,19,400);
S=vector(#R,i,(-1)^(19+i)*R[i]);
print("S_n=(-1)^n R_n all same sign? ",vecmax(vector(#S,i,sign(S[i])!=sign(S[1])))==0);
print("log10|R_n| at n=20,50,100,200,400: ",vector(5,i,my(k=[1,31,81,181,381][i]);round(log(abs(1.0*R[k]))/log(10))));
G=0;
{for(i=1,#R,G=gcd(G,R[i]));}
print("gcd of all R_n digits: ",#Str(G));
DIF=S;
{for(m=1,40, DIF=vector(#DIF-1,i,DIF[i+1]-DIF[i]); my(z=1); for(i=1,#DIF,if(DIF[i]!=0,z=0;break)); if(z,print("  (-1)^n R_n is a POLYNOMIAL in n of degree ",m-1); break); if(m<=20,print("   diff order ",m," max digits ",vecmax(vector(#DIF,i,#Str(abs(DIF[i])))))));}
print("--- ratio test: is R_n hypergeometric? ---");
{for(i=1,6, my(q=R[i+1]/R[i]); print("  R_",20+i,"/R_",19+i," numdig=",#Str(abs(numerator(q)))," dendig=",#Str(abs(denominator(q)))," approx=",1.0*q));}
quit;
