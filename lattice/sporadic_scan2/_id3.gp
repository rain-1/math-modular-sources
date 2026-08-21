default(parisizemax, 6000000000);
default(realprecision, 560);
x = eval(externstr("cat _lim136b.txt")[1]);
print("digits of x read: ", precision(x));
G=znstar(8,1); mf=mfinit([8,3,[G,znconreylog(G,3)]],0); f=mfeigenbasis(mf)[1]; lf=lfunmf(mf,f);
L1=lfun(lf,1); L2=lfun(lf,2);
print("lindep[1,x,L(f8,2)]      = ", lindep([1,x,L2]));
print("lindep[1,x,L(f8,1),L(f8,2)] = ", lindep([1,x,L1,L2]));
print("lindep[1,x,zeta(3)]      = ", lindep([1,x,zeta(3)]));
print("lindep[1,x,zeta(2)]      = ", lindep([1,x,zeta(2)]));
print("lindep[1,x,Pi^2]         = ", lindep([1,x,Pi^2]));
print("lindep[1,x,lfun(-8,2)]   = ", lindep([1,x,lfun(-8,2)]));
print("lindep[1,x,lfun(8,2)]    = ", lindep([1,x,lfun(8,2)]));
print("lindep[1,x,Catalan]      = ", lindep([1,x,Catalan]));
print("lindep[1,x,log(1+sqrt(2))^2] = ", lindep([1,x,log(1+sqrt(2))^2]));
print("lindep[1,x,Pi*log(1+sqrt(2))] = ", lindep([1,x,Pi*log(1+sqrt(2))]));
print("lindep[1,x,L2,zeta(3)]   = ", lindep([1,x,L2,zeta(3)]));
quit
