\\ 07_ctrl.gp -- CONTROL: the known Kohnen form h in S^+_{5/2}(Gamma_0(28)) must satisfy
\\ h*theta^3*Delta^3 in M_40(Gamma_0(28)).  Validates the multiplier bookkeeping.
default(parisize, 4000000000);
NC = 300; PREC = NC+90;
S = mfinit([28,5/2],1); K = mfkohnenbasis(S); h = mflinear(S, K[,1]);
hc = mfcoefs(h, NC+40);
print("h coefs 0..40: ", vector(41,i,hc[i]));
TH = Ser(vector(PREC, i, my(n=i-1); if(n==0,1, if(issquare(n),2,0))), 'q, PREC);
DE = 'q*eta('q+O('q^PREC))^24 + O('q^PREC);
HS = Ser(vector(NC+41, i, hc[i]), 'q, NC+41);
{ for(r=0,3,
   my(w=4+12*r, T = HS*TH^3*DE^r + O('q^(NC+1)), M, B, nb, A, Y, V);
   M = mfinit([28,w],4); B = mfcoefs(M,NC); nb = matsize(B)[2];
   A = matrix(NC+1, nb); for(i=1,NC+1, for(j=1,nb, A[i,j]=B[i,j]));
   Y = vector(NC+1, i, polcoeff(T,i-1));
   V = matinverseimage(A, Y~);
   print("r=",r," weight=",w," dim=",nb," -> ", if(#V==0,"NO","YES: h*theta^3*Delta^r IS in M_w")));
}
quit;
