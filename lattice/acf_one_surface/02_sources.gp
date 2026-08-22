/* ============================================================
   02_sources.gp -- the three sources Phi_X = F_X * D t_X on X_0(6)
   ============================================================ */
read("lib.gp");

E1 = etaprod(1); E2 = etaprod(2); E3 = etaprod(3); E6 = etaprod(6);
E4 = etaprod(4); E12 = etaprod(12);
tC = 'q * E1^4 * E6^8 / (E2^8 * E3^4);
FC = E2^6 * E3 / (E1^3 * E6^2);
tFz = 'q * E1^5 * E3 * E4^5 * E6^2 * E12 / E2^14;
FFz = E2^15 * E3^2 * E12^2 / (E1^6 * E4^6 * E6^5);

tA  = tC/(1-tC);      FA  = (1-tC)*FC;
tFp = tC/(1-9*tC);    FFp = (1-9*tC)*FC;

Es = Ess();  Et = Ett();

PhiC  = FC  * Dop(tC);
PhiA  = FA  * Dop(tA);
PhiFp = FFp * Dop(tFp);
PhiF  = FFz * Dop(tFz);

print("=== A. sources as Eisenstein combinations (level 6, weight 3, chi_{-3}) ===");
print("Phi_C  - (1-8V2) Ess          : ", PhiC  - (Es - 8*Vop(2,Es)));
print("Phi_F' - (1+ V2) Ess          : ", PhiFp - (Es +   Vop(2,Es)));
print("Phi_A  - (1- V2) Ett          : ", PhiA  - (Et -   Vop(2,Et)));
print("Phi_F  - (1-7V2-8V4) Ess      : ", PhiF  - (Es - 7*Vop(2,Es) - 8*Vop(4,Es)));
print("Phi_F  + Phi_F'(tau+1/2)      : ", PhiF  + twist(PhiFp));
print();
print("=== B. the multiplicative (cusp-move) form ===");
print("Phi_A  - Phi_C/(1-t_C)        : ", PhiA  - PhiC/(1-tC));
print("Phi_F' - Phi_C/(1-9 t_C)      : ", PhiFp - PhiC/(1-9*tC));
print();
print("=== C. modular-function identities behind the Eisenstein bookkeeping ===");
print("(1-9 t_C)*(1+V2)Ess - (1-8V2)Ess : ", (1-9*tC)*(Es+Vop(2,Es)) - (Es-8*Vop(2,Es)));
print("(1- t_C)*(1-V2)Ett  - (1-8V2)Ess : ", (1-tC)*(Et-Vop(2,Et)) - (Es-8*Vop(2,Es)));
print();
print("=== D. Eichler integrals Theta = D^{-2} Phi ===");
Ecal  = Dinv2(Es);      \\ D^{-2} S
EcalT = Dinv2(Et);
ThC  = Dinv2(PhiC);  ThFp = Dinv2(PhiFp);  ThA = Dinv2(PhiA); ThF = Dinv2(PhiF);
Psi0 = Es - 4*Vop(2,Es);          ThPsi = Dinv2(Psi0);
print("Th_C  - (1-2V2)Ecal            : ", ThC  - (Ecal - 2*Vop(2,Ecal)));
print("Th_F' - (1+V2/4)Ecal           : ", ThFp - (Ecal + Vop(2,Ecal)/4));
print("Th_A  - (1-V2/4)EcalT          : ", ThA  - (EcalT - Vop(2,EcalT)/4));
print("Th_Psi0 - (1-V2)Ecal           : ", ThPsi - (Ecal - Vop(2,Ecal)));
print("Th_F' + (5/4)Th_C - (9/4)Th_Psi0 : ", ThFp + 5*ThC/4 - 9*ThPsi/4);
print();
print("=== E. Mellin factors P_X(Y) at Y = 1/4 (Y = V_2/4 acting on constants) ===");
print("P_C(Y)=1-8Y  -> P_C(1/4)  = ", 1-8/4);
print("P_F'(Y)=1+Y  -> P_F'(1/4) = ", 1+1/4);
print("P_F(Y)=(1+Y)(1-8Y) -> P_F(1/4) = ", (1+1/4)*(1-8/4));
print("P_A(Y)=1-Y   -> P_A(1/4)  = ", 1-1/4);
quit;
