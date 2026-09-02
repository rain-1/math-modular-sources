\\ second attempt at the three unidentified K_-: allow both Deligne periods (omega^+ ~ L(f,3)/Pi^3, omega^- ~ L(f,2)/Pi^2)
default(realprecision, 70);
KM = [0.20376662415349320849210895209670200943480458841140584624948973756724, 0.37960819633459956111576090173690623070838904894105304002130538432955, 0.43491056992916815669238152561817236709908148635750726222805382955345];
NAMES = ["N10C25B6  Phi=f5-4f5(2t)", "N12C1B34  Phi=f6-4f6(2t)", "N18C1B14  Phi=f9-4f9(2t)"];
LEV = [5, 6, 9];
for(i=1,3, my(M=LEV[i], mf=mfinit([M,4],0), f=mfeigenbasis(mf)[1], L2=lfun(lfunmf(mf,f),2), L3=lfun(lfunmf(mf,f),3), K=KM[i], KS=K*sqrt(Pi)); print("-- ", NAMES[i]); print("     lindep[KS, L2, L3/Pi]         = ", lindep([KS,L2,L3/Pi],45)); print("     lindep[KS, L2, L3/Pi, L2*Pi]  = ", lindep([KS,L2,L3/Pi,L2*Pi],40)); print("     KS/L2 = ", KS/L2, "   algdep8: ", algdep(KS/L2,8)));
quit;
