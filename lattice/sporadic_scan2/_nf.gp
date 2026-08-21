default(realprecision, 40);
G=znstar(8,1);
mf=mfinit([8,3,[G,znconreylog(G,3)]],0);
EB=mfeigenbasis(mf);
print("newdim=",mfdim(mf)," #EB=",#EB);
for(e=1,#EB, print("f",e,"=",mfcoefs(EB[e],16)));
print("eta_1^2 eta_2 eta_4 eta_8^2 :", mfcoefs(mffrometaquo([1,2;2,1;4,1;8,2]),16));
print("eta_2^3 eta_4^3/... test: ", mfcoefs(mffrometaquo([2,3;4,3]),16));
lf=lfunmf(mf,EB[1]); print("L(f,1)=",lfun(lf,1)); print("L(f,2)=",lfun(lf,2)); print("L(f,3)=",lfun(lf,3));
quit
