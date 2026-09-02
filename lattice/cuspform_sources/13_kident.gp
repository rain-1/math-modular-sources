\\ Identification attempts for K_- (15-17 reliable digits from 11_kminus.gp).
default(realprecision, 30);
KM = [0.20376662415349320849, 0.28757331148922160500, 0.37960819633459955517, 0.49809158640016257987, 0.43491056992916815669];
NAMES = ["N10C25B6 f5-4f5(2t)", "N12C9B10 f6-4f6(2t)", "N12C1B34 f6-4f6(2t)", "N18C1B14 f6-9f6(3t)", "N18C1B14 f9-4f9(2t)"];
for(i=1,5, my(K=KM[i]); print("-- ", NAMES[i], "   K_- = ", K); print("     K_-/sqrt(Pi)   = ", K/sqrt(Pi), "     algdep deg4 : ", algdep(K/sqrt(Pi),4)); print("     K_-*sqrt(Pi)   = ", K*sqrt(Pi), "     algdep deg4 : ", algdep(K*sqrt(Pi),4)); print("     K_-*Pi^(3/2)   = ", K*Pi^(3/2), "     algdep deg4 : ", algdep(K*Pi^(3/2),4)); print("     K_-/Pi^(3/2)   = ", K/Pi^(3/2), "     algdep deg4 : ", algdep(K/Pi^(3/2),4)); print("     coefficient of (1-lam2 x)^(1/2) in xi*A-B = -2 sqrt(Pi) K_- = ", -2*sqrt(Pi)*K));
quit;
