default(realprecision, 110);
f = mffrometaquo([2,3;6,3]);
mf = mfinit(f, 0);
lf = lfunmf(mf, f);
L2 = lfun(lf, 2); L3 = lfun(lf,3); L1=lfun(lf,1);
print("L(f,1)=",L1); print("L(f,2)=", L2); print("L(f,3)=", L3);
x = -0.359753349488044326868;
print("lindep 1,x,L2 = ", lindep([1,x,L2]));
