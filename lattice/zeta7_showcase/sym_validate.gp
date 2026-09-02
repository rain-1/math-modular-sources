default(parisize,"4G");
read("/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/analyze_lib.gp");
PV=[n^3, -(34*n^3-51*n^2+27*n-5), (n-1)^3];
v=pol2vec(PV,2,3);
print("apery v = ",v);
{for(A=-6,6, my(w=reflect(v,2,3,A)); if(w==v,print("  A=",A," eps=+1")); if(w==-v,print("  A=",A," eps=-1")));}
print("--- now the mod-p space test on apery ---");
PP=2^61-1;
refmat(A,D)=matrix(D+1,D+1,mm,kk,my(m=mm-1,k=kk-1);if(k>=m,binomial(k,m)*(-1)^m*Mod(A,PP)^(k-m),0));
refvec(v,r,D,TT)=my(w=vector((r+1)*(D+1)));for(j=0,r,my(blk=TT*vector(D+1,i,v[j*(D+1)+i])~);for(m=0,D,w[(r-j)*(D+1)+m+1]=blk[m+1]));w;
Bp=matrix(1,12,a,b,Mod(v[b],PP));
{for(twoA=-12,12, my(A=twoA/2,TT=refmat(A,3),W=matrix(1,12,a,b,0)); my(rv=refvec(Bp[1,],2,3,TT)); for(b=1,12,W[1,b]=rv[b]); if(matrank(matconcat([Bp;W]))==1, print("  space preserved at A=",A)));}
quit;
