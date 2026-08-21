default(parisize,2000000000); default(realprecision,60);
PRECSET=900; read("lattice/rigidity/lib.gp");
tF = etaq([[1,5],[3,1],[4,5],[6,2],[12,1],[2,-14]]); FF = etaq([[2,15],[3,2],[12,2],[1,-6],[4,-6],[6,-5]]);
PhiF = FF*Dq(tF); ThF = Dqinv(PhiF,2);
ev(f,tau)={my(Q=exp(2*Pi*I*tau),s=0,P=1); for(n=0,PRECSET-1, s+=polcoeff(f,n)*P; P*=Q); s};
gam(t)=(8*t-3)/(12*t-4);
Pp(T)=(6*T-2)*ev(ThF,gam(T))-ev(ThF,T);
xinf = lfun(-3,2)/2;
Pth(T) = 5/4*xinf*(6*T-3) - I*Pi^2/(6*sqrt(3))*(6*T-1);
tst(T)=print(T,"   err=",abs(Pp(T)-Pth(T)));
tst(0.37+0.31*I); tst(0.30+0.26*I); tst(0.42+0.35*I); tst(0.36+0.40*I); tst(0.55+0.29*I); tst(0.25+0.5*I); tst(-0.1+0.45*I);
print("xi_inf(C)=L(chi-3,2)/2 = ",xinf);
quit
