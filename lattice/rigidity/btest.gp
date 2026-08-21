PS=60; read("cover.gp"); read("relfind.gp");
\\ basis of Q(u,v) over Q: u^i v^j, j<=1 (v^2 reduces via modular equation)
tryin(g,nm,di)={my(LL=List(),TT=List(),K);
 for(i=0,di,for(j=0,1, listput(LL,u^i*v^j); listput(TT,[0,i,j]); listput(LL,g*u^i*v^j); listput(TT,[1,i,j])));
 K=relker(Vec(LL),50); print(nm," (deg ",di,") kerdim ",#K);
 if(#K>0, for(c=1,min(#K,2), my(s=K[,c]); s=s/content(s); print("   ",sum(i=1,#TT, s[i]*'T^TT[i][1]*'U^TT[i][2]*'V^TT[i][3]))));}
tryin(tB,"tB",6);

w3=Vd(tC,3);
tryin(w3,"tC(3tau)",6);
\\ relation between tB and tF
LL=List();TT=List();
for(i=0,4,for(j=0,4, listput(LL,tB^i*tF^j); listput(TT,[i,j])));
K=relker(Vec(LL),50); print("tB vs tF kerdim ",#K);
if(#K>0, my(s=K[,1]); s=s/content(s); print("   ",sum(i=1,#TT, s[i]*'B^TT[i][1]*'X^TT[i][2])));
\\ relation between tB and tC
LL=List();TT=List();
for(i=0,4,for(j=0,4, listput(LL,tB^i*tC^j); listput(TT,[i,j])));
K=relker(Vec(LL),50); print("tB vs tC kerdim ",#K);
if(#K>0, my(s=K[,1]); s=s/content(s); print("   ",sum(i=1,#TT, s[i]*'B^TT[i][1]*'U^TT[i][2])));
