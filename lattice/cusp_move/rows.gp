/* rows.gp -- the corpus of second-order Apery-like rows.
   Each entry: [P(nv), Q(nv), "name", "provenance"].
   Convention (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1}, u_0=1, u_1=P(0).  */

{
corpus = [
 /* --- Zagier's six (Q = d n^2, r1=r2=0) --- */
 [7*nv^2+7*nv+2,    -8*nv^2,               "ZagierA",  "(7,2,-8)   zeta(2)/4"],
 [9*nv^2+9*nv+3,    27*nv^2,               "ZagierB",  "(9,3,27)   no arch limit"],
 [10*nv^2+10*nv+3,  9*nv^2,                "ZagierC",  "(10,3,9)   L(2,chi-3)/2"],
 [11*nv^2+11*nv+3,  -nv^2,                 "ZagierD",  "(11,3,-1)  zeta(2)/5"],
 [12*nv^2+12*nv+4,  32*nv^2,               "ZagierE",  "(12,4,32)  G/2 Catalan"],
 [17*nv^2+17*nv+6,  72*nv^2,               "ZagierF",  "(17,6,72)  5L(2,chi-3)/8"],
 /* --- root rows, Q = c(2n-1)^2, r1=r2=1/2 --- */
 [136*nv^2+68*nv+10,  4*(2*nv-1)^2,        "sqrtApery","Beukers (136,10,4) = sqrt(Apery)"],
 [24*nv^2+12*nv+2,    4*(2*nv-1)^2,        "sqrtT",    "sqrt of T=(12,4,16)"],
 [20*nv^2+10*nv+2,    16*(2*nv-1)^2,       "sqrtDomb", "sqrt of Domb (10,4,64)"],
 [72*nv^2+36*nv+6,    -108*(2*nv-1)^2,     "sqrtAZ93", "sqrt of AZ(9,3,-27)"],
 [88*nv^2+44*nv+10,   500*(2*nv-1)^2,      "sqrtAZ11", "sqrt of AZ(11,5,125)"],
 [56*nv^2+28*nv+6,    324*(2*nv-1)^2,      "sqrtAZ73", "sqrt of AZ(7,3,81)"],
 /* --- Cooper root rows (Q with distinct roots) --- */
 [26*nv^2+13*nv+2,    -3*(3*nv-1)*(3*nv-2),"sqrts7",   "sqrt of Cooper s7"],
 [24*nv^2+12*nv+2,    -4*(8*nv-3)*(8*nv-5),"sqrts10",  "sqrt of Cooper s10"],
 [56*nv^2+28*nv+6,    12*(8*nv-3)*(8*nv-5),"sqrts18",  "sqrt of Cooper s18"],
 /* --- the Gamma_0(5)+5 non-congruence row --- */
 [88*nv^2+44*nv+6,    -4*(4*nv-1)*(4*nv-3),"level5F",  "Gamma_0(5)+5 non-congruence"],
 /* --- the two new Herfurtner rows --- */
 [117*nv^2+78*nv+21,  441*(3*nv-1)^2,      "herf30",   "#30 I1 I7 II II"],
 [72*nv^2+36*nv+6,    108*(4*nv-1)*(4*nv-3),"herf45",  "#45 I3 III III III"]
];
}
