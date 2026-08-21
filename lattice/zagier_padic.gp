\\ Zagier sporadic second-order families: p-adic slopes and cross-row coincidences.
default(realprecision,60);
fam=["A","B","C","D","E","F"];
par=[[7,2,-8],[9,3,27],[10,3,9],[11,3,-1],[12,4,32],[17,6,72]];
N=260;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;
            B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
R=vector(6,i,rows(par[i][1],par[i][2],par[i][3]));
print("family  c   limit b_n/a_n   max v_p(a_n), p=2,3,5,7");
{for(i=1,6,
 my(A=R[i][1],B=R[i][2],c=par[i][3],lim=B[N+1]/A[N+1]);
 my(mv=vector(4,j,my(p=[2,3,5,7][j]);vecmax(vector(N,n,valuation(A[n+1],p)))));
 print(fam[i],"  ",c,"  ",lim,"  ",mv));}
print("\ncross-row v_p(a^X_n b^Y_n - a^Y_n b^X_n) at n=130 and 260, p=2 | p=3");
{for(i=1,6,for(j=i+1,6,
 my(A1=R[i][1],B1=R[i][2],A2=R[j][1],B2=R[j][2]);
 my(d1=A1[131]*B2[131]-A2[131]*B1[131], d2=A1[261]*B2[261]-A2[261]*B1[261]);
 print(fam[i],fam[j],"  p=2: ",valuation(d1,2)," ",valuation(d2,2),"   p=3: ",valuation(d1,3)," ",valuation(d2,3))));}
\q
