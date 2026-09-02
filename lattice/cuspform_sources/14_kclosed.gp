\\ Closed form of K_-: search for  K_- sqrt(Pi) = a L(f,2) + b L(f,3)  over Q and Q(sqrt d).
default(realprecision, 70);
KM = [0.20376662415349320849210895209670200943480458841140584624948973756724, 0.28757331148922160500069983443231806149556558391059014893586330391409, 0.37960819633459956111576090173690623070838904894105304002130538432955, 0.49809158640016257987401106018289485845958492345232739269018567197489, 0.43491056992916815669238152561817236709908148635750726222805382955345];
NAMES = ["N10C25B6  Phi=f5-4f5(2t)", "N12C9B10  Phi=f6-4f6(2t)", "N12C1B34  Phi=f6-4f6(2t)", "N18C1B14  Phi=f6-9f6(3t)", "N18C1B14  Phi=f9-4f9(2t)"];
LEV = [5, 6, 6, 6, 9];
for(i=1,5, my(M=LEV[i], mf=mfinit([M,4],0), f=mfeigenbasis(mf)[1], L1=lfun(lfunmf(mf,f),1), L2=lfun(lfunmf(mf,f),2), L3=lfun(lfunmf(mf,f),3), K=KM[i], KS=K*sqrt(Pi)); print("-- ", NAMES[i], "   newform level ", M); print("     K_- sqrt(Pi)              = ", KS); print("     K_- sqrt(Pi) / L(f,2)     = ", KS/L2, "    algdep4 : ", algdep(KS/L2,4)); print("     K_- sqrt(Pi) / L(f,3)     = ", KS/L3, "    algdep4 : ", algdep(KS/L3,4)); print("     lindep[KS, L2, L3]        = ", lindep([KS,L2,L3],45)); print("     lindep[KS, L2, L3, Pi^3]  = ", lindep([KS,L2,L3,Pi^3],40)));
quit;
