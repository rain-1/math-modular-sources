default(realprecision, 210);
LK(D, s) = { my(q=abs(D), t=0); for(j=1,q-1, t += kronecker(D,j)*zetahurwitz(s, j/q)); t/q^s; }
mf7 = mfinit([7,3,-7],0); f7 = mfeigenbasis(mf7)[1]; L7 = lfunmf(mf7,f7);
G12 = znstar(12,1); mf12 = mfinit([12,3,[G12,[0,1]]],0); f12 = mfeigenbasis(mf12)[1]; L12 = lfunmf(mf12,f12);
mf27 = mfinit([27,3,[znstar(27,1),[9]]],0); E27 = mfeigenbasis(mf27);
print("dim27=",#E27); for(i=1,#E27, print("  f27_",i," ",mfcoefs(E27[i],10)));
write("lvals.txt","Lf7_2 ", lfun(L7,2));
write("lvals.txt","Lf7_1 ", lfun(L7,1));
write("lvals.txt","Lf12_2 ", lfun(L12,2));
write("lvals.txt","Lf12_1 ", lfun(L12,1));
write("lvals.txt","L7c ", LK(-7,2));
write("lvals.txt","L3c ", LK(-3,2));
write("lvals.txt","Cat ", Catalan);
write("lvals.txt","Pi ", Pi);
write("lvals.txt","z3 ", zeta(3));
