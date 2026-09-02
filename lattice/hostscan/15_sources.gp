/* 15_sources.gp -- the companion source Phi = F * theta_q(x) for each of the 12
   CDT-shape Fricke rows: is it holomorphic (Eisenstein) on Gamma_0(N)?  If so,
   express it in the Eisenstein basis of M_4(Gamma_0(N)) and record the level. */
default(parisizemax, 8000000000);
NQ = 60;
useries(dv, r, nq) = q*prod(t=1, #dv, eta(q^dv[t] + O(q^nq))^r[t]);
Fseries(dv, r, nq) = 1 - sum(t=1, #dv, r[t]*dv[t]*sum(n=1, (nq-1)\dv[t], sigma(n)*q^(dv[t]*n))) + O(q^nq);
thq(s) = my(n=serprec(s,q)); sum(k=0,n-1, k*polcoeff(s,k)*q^k) + O(q^n);
doit(nn,dv,r,C,B,nm) = my(us=useries(dv,r,NQ), Fs=Fseries(dv,r,NQ)); my(xs=us/(1+B*us+C*us^2)); my(Ph=Fs*thq(xs)); my(v=vector(NQ-2,i,polcoeff(Ph,i-1))); print("== ",nm); print("   Phi = ", vector(9,i,polcoeff(Ph,i-1))); my(mf=mfinit([nn,4],4)); my(cf=mftobasis(mf,v,1)); print("   dim M_4(Gamma_0(",nn,")) = ", mfdim(mf), "   mftobasis residual-flag coords: ", cf); my(mfe=mfinit([nn,4],3)); print("   dim M_4^Eis = ", mfdim(mfe), "   in Eis? ", mftobasis(mfe,v,1));
doit(5,[1,5],[-6,6],125,22,"N=5 C=125 B=22  (zeta(2)/10)");
doit(6,[1,2,3,6],[-4,-4,4,4],81,14,"N=6 C=81 B=14  (zeta(2)/8)");
doit(6,[1,2,3,6],[-5,1,-1,5],72,17,"N=6 C=72 B=17  (Apery zeta(3)/6)");
doit(6,[1,2,3,6],[-6,6,-6,6],64,20,"N=6 C=64 B=20  (L(2,chi-3)/4)");
doit(7,[1,7],[-4,4],49,13,"N=7 C=49 B=13  (Cooper s7, zeta(2)/7)");
doit(8,[1,2,4,8],[-4,2,-2,4],32,12,"N=8 C=32 B=12  (AZ eps, 7zeta(3)/32)");
doit(8,[1,2,4,8],[-8,16,-16,8],16,24,"N=8 C=16 B=24  (G/4)");
doit(9,[1,3,9],[-3,0,3],27,9,"N=9 C=27 B=9  (AZ zeta, L(3,chi-3)/3)");
doit(10,[1,2,5,10],[-2,-2,2,2],25,6,"N=10 C=25 B=6  (Cooper s10)");
doit(12,[1,2,3,4,6,12],[-4,4,4,-4,-4,4],9,10,"N=12 C=9 B=10  (Domb)");
doit(12,[1,2,3,4,6,12],[-12,24,12,-12,-24,12],1,34,"N=12 C=1 B=34  (5L(2,chi-3)/16)");
doit(18,[1,2,3,6,9,18],[-6,6,12,-12,-6,6],1,14,"N=18 C=1 B=14  (Cooper s18)");
quit;
