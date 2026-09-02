default(parisize,"8G");
read("lib.gp");
NT = 800;
IDX = [8, 4];
{ for(u=1,2,
  my(i=IDX[u], R=ROWS[i], r=R[2], AB=genrow(R,NT), A=AB[1], B=AB[2], dn=vector(NT+1), NN=vector(NT+1));
  dn[1]=1; for(n=1,NT, dn[n+1]=lcm(dn[n],n));
  for(n=0,NT, NN[n+1] = dn[n+1]^r*B[n+1]);
  print("\n===== row ", R[1], " r=", r, " : N_n = d_n^r b_n ; N_0..N_8 = ", vector(9,j,NN[j]), " =====");
  for(pi=1,4,
    my(p=[5,7,11,13][pi], okA=0, totA=0, okN=0, totN=0, nmax=floor(NT/p)-1);
    for(n=1,nmax, for(m=0,p-1,
      totA++; if((A[n*p+m+1] - A[n+1]*A[m+1])%p==0, okA++);
      totN++; if((NN[n*p+m+1] - NN[n+1]*NN[m+1])%p==0, okN++)));
    print("  p=", p, " Lucas a_{np+m}=a_n a_m mod p : ", okA, "/", totA, "   Lucas N_{np+m}=N_n N_m mod p : ", okN, "/", totN);
    my(rt=List(), vA=List(), vN=List());
    for(n=1,min(24,nmax), 
      if(NN[n+1]%p!=0, listput(rt, lift(Mod(NN[n*p+1],p)/Mod(NN[n+1],p))), listput(rt,-1));
      listput(vA, if(A[n*p+1]==A[n+1], 99, valuation(A[n*p+1]-A[n+1],p)));
      listput(vN, if(NN[n*p+1]==NN[n+1], 99, valuation(NN[n*p+1]-NN[n+1],p))));
    print("      N_{np}/N_n mod p, n=1..24 : ", Vec(rt));
    print("      v_p(a_{np} - a_n),  n=1..24 : ", Vec(vA));
    print("      v_p(N_{np} - N_n),  n=1..24 : ", Vec(vN))));
}
quit;
