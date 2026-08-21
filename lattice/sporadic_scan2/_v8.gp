default(realprecision, 70);
G=znstar(8,1);
mf=mfinit([8,3,[G,znconreylog(G,3)]],0);
f=mfeigenbasis(mf)[1];
lf=lfunmf(mf,f);
print("L(f,2)=", lfun(lf,2));
print("f = ", mfcoefs(f,12));
print("eta_1^2 eta_2 eta_4 eta_8^2 = ", mfcoefs(mffrometaquo([1,2;2,1;4,1;8,2]),12));
\\ the weight-1 form: theta series?
mf1=mfinit([8,1,[G,znconreylog(G,3)]],4);
print("dim M_1(8,chi_-8)=", mfdim(mf1));
print("basis coefs:", mfcoefs(mf1,16));
quit
