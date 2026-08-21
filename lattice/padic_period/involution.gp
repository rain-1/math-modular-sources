\\ The involution gamma = [8,-3;12,-4] (Atkin-Lehner W_4 on X_0(12)) realising w.
default(parisize,2000000000); default(realprecision,60);
PRECSET=600; read("lattice/rigidity/lib.gp");
tC = etaq([[1,4],[6,8],[2,-8],[3,-4]]); FC = etaq([[2,6],[3,1],[1,-3],[6,-2]]);
tF = etaq([[1,5],[3,1],[4,5],[6,2],[12,1],[2,-14]]); FF = etaq([[2,15],[3,2],[12,2],[1,-6],[4,-6],[6,-5]]);
PhiF = FF*Dq(tF);
ev(f,tau)={my(Q=exp(2*Pi*I*tau),s=0,P=1); for(n=0,PRECSET-1, s+=polcoeff(f,n)*P; P*=Q); s};
gam(t)=(8*t-3)/(12*t-4);  \\ normalised: [4,-3/2;6,-2] in SL_2(R), j = 6*tau-2
sig(t)=(t-1)/(9*t-1);
one(T)={my(g=gam(T),j=6*T-2,u,v);u=ev(tC,T);v=ev(tC,2*T);
 print("tau=",T,"  Im(gam tau)=",imag(g));
 print("  tC(g)-ubar   = ",abs(ev(tC,g)-sig(v)));
 print("  tC(2g)-vbar  = ",abs(ev(tC,2*g)-sig(u)));
 print("  tF(g)-tF     = ",abs(ev(tF,g)-ev(tF,T)));
 print("  FF|_1 g / FF = ",ev(FF,g)/(j*ev(FF,T)));
 print("  PhiF|_3 g/PhiF = ",ev(PhiF,g)/j^3/ev(PhiF,T));};
one(0.37+0.31*I); one(0.30+0.26*I); one(0.42+0.35*I);
quit
