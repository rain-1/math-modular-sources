\\ 05_find.gp -- Is  F = sum_d c(d) q^d  (c from round2's twisted CM traces, plus an unknown
\\ polar coefficient at d=-3 and unknowns at the DEGEN indices) a weakly holomorphic form of
\\ weight 5/2 on Gamma_0(28)?  Test: F*theta^3*Delta^3 must lie in M_40(Gamma_0(28)).
default(parisize, 4000000000);
NC = 600;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(2500);
KNOWN = vector(2500);
DEG = List();
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2]);
    if(d>2500, next);
    if(type(v)=="t_STR", listput(DEG,d), CD[d]=v; KNOWN[d]=1));
}
print("entries=", #DAT, "  DEGEN=", Vec(DEG));
\\ q-series
T3 = Ser(vector(NC+40, i, my(n=i-1); sum(a=-40,40, if(a^2==n,1,0))), 'q, NC+40)^3;
D3 = (Ser(vector(NC+40,i,polcoeff(eta('q + O('q^(NC+40)))^24,i-1)),'q,NC+40))^3;
P = T3*D3;   \\ starts q^3
print("P valuation = ", valuation(P,'q));
\\ unknown indices: -3 and the DEGEN d <= NC
UN = List(); listput(UN,-3);
{ for(i=1,#DEG, if(DEG[i]<=NC, listput(UN,DEG[i]))); }
UN = Vec(UN);
print("unknowns at d = ", UN);
\\ known part K
KS = sum(d=1, NC, if(KNOWN[d], CD[d], 0)*'q^d) + O('q^(NC+1));
RHS = KS*P;   \\ series
\\ basis of M_40(28)
gettime();
M = mfinit([28,40],4);
B = mfcoefs(M, NC);
print("dim=", matsize(B), " time=", gettime());
\\ build matrix: columns = basis (159) then -q^{d_j}*P for each unknown
nb = matsize(B)[2];
nu = #UN;
A = matrix(NC+1, nb+nu);
{ for(i=1,NC+1, for(j=1,nb, A[i,j]=B[i,j])); }
{ for(j=1,nu, my(S = 'q^UN[j]*P + O('q^(NC+1)));
    for(i=1,NC+1, A[i,nb+j] = -polcoeff(S, i-1))); }
Y = vector(NC+1, i, polcoeff(RHS, i-1));
print("matrix size ", matsize(A), "  rank=", matrank(A));
V = matinverseimage(A, Y~);
if(#V==0, print("NO SOLUTION -- F is NOT in the space"); quit);
RES = A*V - Y~; print("residual over all ", NC+1, " rows: ", if(RES==0, "ZERO -- F IS MODULAR", "NONZERO"));
if(RES!=0, print("nonzero residual rows: ", select(x->x, vector(NC+1,i,if(RES[i]!=0,i-1,0)), 1)));
print("polar coefficient c(-3) = ", V[nb+1]);
{ for(j=2,nu, print("c(", UN[j], ") = ", V[nb+j])); }
KER = matker(A);
print("kernel dim = ", matsize(KER)[2]);
write("05_sol.txt", V);
quit;
