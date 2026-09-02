default(parisize,"20G");
read("/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/analyze_lib.gp");
pencil(fn,r,D)=my(B,k,CH,G,best=[-1,-1],bestv=0);read(concat(BASE,fn));B=BVL;k=#B;CH=vector(k,i,chiof(B[i],r,D));print("---- (",r,",",D,") dim ",k," ----");G=CH[1];for(i=2,k,G=gcd(G,CH[i]));G=G/content(G);print("  gcd of chi over basis: degree ",poldegree(G));if(poldegree(G)>0,print("  gcd factored: ",factor(G)));for(a=0,r,for(b=0,(r-a)\2,my(md=(L+1)^a*(L^2-14*L+1)^b,MM,kk);if(a+2*b<=r,MM=matrix(a+2*b,k,x,y,polcoeff(CH[y]%md,x-1));kk=matker(MM*1.0*0+MM);if(#kk>0,if(a+2*b>best[1]+2*best[2],best=[a,b];bestv=kk)))));print("  max (a,b) with (L+1)^a*(L^2-14L+1)^b | chi : ",best,"  freedom dim ",if(type(bestv)=="t_MAT",matsize(bestv)[2],0));return([B,CH,best,bestv]);
gettime();
R1=pencil("rec_23_10_raw.gp",23,10);
print("time ",gettime());
R2=pencil("rec_17_15_raw.gp",17,15);
print("time ",gettime());
quit;
