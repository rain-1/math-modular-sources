mf=mfinit([5,4,0],3);
print("params=",mfparams(mf));
B=mfbasis(mf);
print("nbasis=",#B);
for(i=1,#B, print(i," ",mfcoefs(B[i],6)));
print("dimvec=",mfdim([5,4,0],3));
print("dimtot=",mfdim(mf));
quit;
