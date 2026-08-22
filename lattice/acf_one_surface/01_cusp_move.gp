/* ============================================================
   01_cusp_move.gp  --  the three Hauptmoduls on X_0(6)
   ============================================================ */
read("lib.gp");

E1 = etaprod(1); E2 = etaprod(2); E3 = etaprod(3); E6 = etaprod(6);
E4 = etaprod(4); E12 = etaprod(12);

/* t_C = eta_1^4 eta_6^8 / (eta_2^8 eta_3^4) = q * E1^4 E6^8 /(E2^8 E3^4) */
tC = 'q * E1^4 * E6^8 / (E2^8 * E3^4);
/* F_C = eta_2^6 eta_3 /(eta_1^3 eta_6^2) */
FC = E2^6 * E3 / (E1^3 * E6^2);

/* Zagier's F, level 12 (CONJ_D_PROOF Table):
   t_F = eta_1^5 eta_3 eta_4^5 eta_6^2 eta_12 / eta_2^14,
   F_F = eta_2^15 eta_3^2 eta_12^2 /(eta_1^6 eta_4^6 eta_6^5)                */
tF = 'q * E1^5 * E3 * E4^5 * E6^2 * E12 / E2^14;
FF = E2^15 * E3^2 * E12^2 / (E1^6 * E4^6 * E6^5);

NR = 240;
uA = zrow(7,2,-8,NR);  uC = zrow(10,3,9,NR);  uFz = zrow(17,6,72,NR);
uFp = vector(NR+1, k, (-1)^(k-1)*uFz[k]);   \\ F' row: u_n^{F'} = (-1)^n u_n^F

print("=== A. rows (first terms) ===");
print("A  : ", vector(8,k,uA[k]));
print("C  : ", vector(8,k,uC[k]));
print("F  : ", vector(8,k,uFz[k]));
print("F' : ", vector(8,k,uFp[k]));

print();
print("=== B. Hauptmodul identities on X_0(6) ===");
tA  = tC/(1-tC);
tFp = tC/(1-9*tC);
print("t_C  = ", tC + O('q^9));
print("t_A  = t_C/(1-t_C)  = ", tA + O('q^9));
print("t_F' = t_C/(1-9t_C) = ", tFp + O('q^9));
print("t_F  = ", tF + O('q^9));
print("check  t_F(tau)  = -t_F'(tau+1/2) :  ", tF + twist(tFp));
print("check  t_F'      = -t_F(tau+1/2)  :  ", tFp + twist(tF));

print();
print("=== C. the weight-one forms ===");
FA  = evalrow(uA , tA , NQ);
FCc = evalrow(uC , tC , NQ);
FFp = evalrow(uFp, tFp, NQ);
FFz = evalrow(uFz, tF , NQ);
vsh(nm,f) = print(nm, " : vanishes to order ", if(f==0, "infinity", valuation(f,'q)), "   (row truncation NR=", NR, ")");
vsh("F_C  from row  vs eta quotient", FCc - FC);
vsh("F_F  from row  vs eta quotient", FFz - FF);
print("F_A  = ", FA + O('q^9));
print("F_F' = ", FFp + O('q^9));
vsh("check F_F'(tau) = F_F(tau+1/2)", FFp - twist(FF));

print();
print("=== D. gauge factors g = F_X / F_C as rational functions of t_C ===");
gA  = FA /FC;
gFp = FFp/FC;
print("g_A  = F_A /F_C = ", gA  + O('q^9));
print("g_F' = F_F'/F_C = ", gFp + O('q^9));
/* express in t_C ; the answer is exactly linear */
gA_t  = texp(gA , tC);
gFp_t = texp(gFp, tC);
print("g_A  as a series in t_C   : ", gA_t  + O('q^6));
print("g_F' as a series in t_C   : ", gFp_t + O('q^6));
vsh("g_A  - (1 - t_C)         ", gA  - (1-tC));
vsh("g_F' - (1 - 9 t_C)       ", gFp - (1-9*tC));
print();
print("=== E. so:  F_A = (1-t_C) F_C ,  F_F' = (1-9t_C) F_C   [EXACT] ===");
quit;
