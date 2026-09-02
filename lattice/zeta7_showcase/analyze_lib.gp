BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
tovec2pol(v,r,D)=vector(r+1,jj,sum(k=0,D,v[(jj-1)*(D+1)+k+1]*n^k));
pol2vec(PV,r,D)=my(w=vector((r+1)*(D+1)));for(jj=0,r,for(k=0,D,w[jj*(D+1)+k+1]=polcoeff(PV[jj+1],k,n)));w;
prim(v)=my(dd=1,g=0);for(i=1,#v,dd=lcm(dd,denominator(v[i])));v=v*dd;for(i=1,#v,g=gcd(g,v[i]));if(g==0,return(v));v/g;
maxdig(v)=my(m=0);for(i=1,#v,m=max(m,#Str(abs(v[i]))));m;
chiof(v,r,D)=sum(jj=0,r,v[jj*(D+1)+D+1]*L^(r-jj));
reflect(v,r,D,A)=my(PV=tovec2pol(v,r,D),QV=vector(r+1,jj,subst(PV[r+2-jj],n,A-n)));pol2vec(QV,r,D);
