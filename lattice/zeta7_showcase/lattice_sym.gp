default(parisize,"20G");
read("/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/analyze_lib.gp");
propq(w,v)=my(i0=0,lam);for(i=1,#v,if(v[i]!=0&&i0==0,i0=i));lam=w[i0]/v[i0];if(w==lam*v,return(lam),return(0));
doshape(fn,r,D)=my(k,B,S,U,RB,hts);read(concat(BASE,fn));B=BVL;k=#B;print("---- shape (",r,",",D,") dim ",k," ----");S=matrixqz(mattranspose(matrix(k,#B[1],a,b,B[a][b])),-1);U=qflll(S);RB=S*U;print("  LLL-reduced basis max digits: ",vector(k,i,maxdig(RB[,i]~)));return(RB);
