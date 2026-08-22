\\ cooper_scores.gp -- Cooper s7,s10,s18: denominator exponent k, kappa_p, sigma_p, lambda_1.
\\ Recurrences from consolidation/SLOPE_CENSUS.md sec.1 (initial data A_1 = 2,4,6).
default(parisizemax,"8G");
default(realprecision,60);
N = 420;

fs10(n) = 2*(2*n+1)*(3*n^2+3*n+1);
gs10(n) =  4*n*(16*n^2-1);
fs7(n)  = (2*n+1)*(13*n^2+13*n+4);
gs7(n)  =  3*n*(9*n^2-1);
fs18(n) = 2*(2*n+1)*(7*n^2+7*n+3);
gs18(n) = -12*n*(16*n^2-1);

{ rowc(idx,N) =
  my(A=vector(N+1), Cc=vector(N+1), a1=[2,4,6][idx], fv, gv);
  A[1]=1; A[2]=a1; Cc[1]=0; Cc[2]=1;
  for(n=1,N-1,
    fv = if(idx==1, fs10(n), idx==2, fs7(n), fs18(n));
    gv = if(idx==1, gs10(n), idx==2, gs7(n), gs18(n));
    A[n+2]=(fv*A[n+1]+gv*A[n])/(n+1)^3;
    Cc[n+2]= fv*Cc[n+1]+gv*n^3*Cc[n];
  );
  [A,Cc];
}
bv(R,n) = R[2][n+1]/(n!)^3;

{ for(idx=1,3,
  my(nm=["s10","s7","s18"][idx], R=rowc(idx,N), d=1, kmax=0, kAmax=0);
  \\ denominator exponents: for A_n (should be 0 = integral) and for B_n
  for(n=1,200, d=lcm(d,n);
    my(kk=0); while(denominator(bv(R,n)*d^kk)!=1 && kk<15, kk++); if(kk>kmax,kmax=kk);
    my(k2=0); while(denominator(R[1][n+1]*d^k2)!=1 && k2<15, k2++); if(k2>kAmax,kAmax=k2));
  print("\n", nm, ":  k(B_n) = ", kmax, "   k(A_n) = ", kAmax);
  print("   log max(|A_n|,|B_n|)/n at n=280,350,420 : ",
    vector(3,j, my(n=[280,350,420][j]); log(max(abs(R[1][n+1]*1.),abs(bv(R,n)*1.)))/n));
  my(xi=bv(R,N)/R[1][N+1]);
  for(j=1,4, my(pp=[2,3,5,7][j],
      v1=valuation(xi-bv(R,200)/R[1][201],pp), v2=valuation(xi-bv(R,400)/R[1][401],pp));
    print("   p=",pp,"  v_p at n=200,400: ",v1,", ",v2,"   slope=",(v2-v1)/200.,
          "   kappa_p=", -valuation(R[1][401],pp)/400.));
); }
quit;
