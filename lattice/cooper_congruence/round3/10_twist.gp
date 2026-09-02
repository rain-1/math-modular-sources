\\ 10_twist.gp -- the s7 vector-valued input is ANTIsymmetric in beta (09_beta.gp).
\\ Canonical scalar avatar:  ct(d) = (beta_d/7) * c(beta_d,d),  independent of the sign of beta.
\\ Test: ct * theta^3 * Delta^3  in  M_40(Gamma_0(28), chi)  for every even chi mod 28.
default(parisize, 8000000000);
NC = 400; PREC = NC+120;
DAT = read("../round2/73_cd_s7_2500.txt");
CD = vector(2600); BU = vector(2600);
{ for(i=1,#DAT, my(e=DAT[i], d=e[1], v=e[2], bs=List(), bt);
    for(b=0,13, if((b^2+3*d)%28==0, listput(bs,b)));
    bs = Vec(bs);
    bt = if(issquare(d), (5*sqrtint(d))%14, bs[1]);
    if(!setsearch(Set(bs), bt), bt = bs[1]);
    BU[d] = bt;
    CD[d] = if(type(v)=="t_STR", 0, v)); }
\\ canonical twist
CT = vector(2600, d, kronecker(BU[d],7)*CD[d]);
print("ct(d) for the first admissible d: ", select(x->x!=0, vector(120,d,CT[d]), 1));
print("ct at d=1,4,8,9,16,25,29,32,36,37: ", vector(10,i,CT[[1,4,8,9,16,25,29,32,36,37][i]]));
TH = Ser(vector(PREC, i, my(n=i-1); if(n==0,1, if(issquare(n),2,0))), 'q, PREC);
DE = 'q*eta('q+O('q^PREC))^24 + O('q^PREC);
S0 = sum(d=1, NC+60, CT[d]*'q^d) + O('q^(NC+61));
G = znstar(28,1);
CYC = G.cyc;
{ run(r, sgn) =
  my(w=4+12*r, PW = TH^3*DE^r, RHS, UN, nu, A, Y, V, M, B, nb, S1);
  S1 = sum(d=1, NC+60, sgn(d)*CT[d]*'q^d) + O('q^(NC+61));
  RHS = S1*PW + O('q^(NC+1));
  UN = [-3, 0];  nu = #UN;
  forvec(ch = vector(#CYC, i, [0, CYC[i]-1]),
    my(CHI = [G, ch~]);
    if(zncharisodd(G, ch~), next);
    M = mfinit([28,w,CHI], 4);
    B = mfcoefs(M, NC); nb = matsize(B)[2];
    if(nb==0, print("  chi=",ch," dim 0"); next);
    A = matrix(NC+1, nb+nu);
    for(i=1,NC+1, for(j=1,nb, A[i,j]=B[i,j]));
    for(j=1,nu, my(S='q^UN[j]*PW + O('q^(NC+1))); for(i=1,NC+1, A[i,nb+j]=-polcoeff(S,i-1)));
    Y = vector(NC+1, i, polcoeff(RHS,i-1));
    V = matinverseimage(A, Y~);
    print("  r=",r," chi=",ch," dim=",nb," -> ", if(#V==0,"NO", concat("YES  c(-3)=",concat(Str(V[nb+1]), concat("  c(0)=",Str(V[nb+2])))))));
}
print("=== plain ct ===");
run(3, d->1);
quit;
