default(realprecision,40);
mf=mfinit([5,3,Mod(2,5)],3);
B=mfbasis(mf);
print("dim=",#B);
print("c0 exact=",vector(#B,i,mfcoefs(B[i],0)[1]));
print("slashI flrat0=",vector(#B,i,mfslashexpansion(mf,B[i],[1,0;0,1],0,0)));
print("slash0 flrat0=",vector(#B,i,mfslashexpansion(mf,B[i],[0,-1;1,0],0,0)));
print("embed=",mfembed(mf,B[1]));
quit;
