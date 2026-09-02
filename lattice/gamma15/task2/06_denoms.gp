/* 06_denoms.gp -- sharp denominator type for BD and B_new in Z[phi].
   Z[phi] = {(u+v sqrt5)/2 : u,v in Z, u = v mod 2};  B_new,n = (u_n+v_n sqrt5)/2,
   u_n = 2*B3_n + 11*B4_n,  v_n = 5*B4_n.                                    */
read("build.txt");
expo(r, d) = my(k=0); while(denominator(r*d^k)>1 && k<15, k++); k;
dn = vector(NA+1); dn[1]=1; for(n=1,NA, dn[n+1]=lcm(dn[n],n));
vD = vector(NA); vN = vector(NA); v3 = vector(NA); v4 = vector(NA); par = vector(NA);
for(n=1, NA, my(d=dn[n+1]); my(u=2*B3[n+1]+11*B4[n+1]); my(v=5*B4[n+1]); vD[n]=expo(BD[n+1],d); v3[n]=expo(B3[n+1],d); v4[n]=expo(B4[n+1],d); vN[n]=max(expo(u,d),expo(v,d)); par[n]=if(vN[n]<=2, ((u*d^2-v*d^2)%2==0), -1));
print("N = ", NA);
print("max_n e(BD)   = ", vecmax(vD), "    max_n e(B3) = ", vecmax(v3), "   max_n e(B4) = ", vecmax(v4), "   max_n e(Bnew u,v) = ", vecmax(vN));
print("all parities u=v mod 2 after d_n^2 ? ", vecmin(par)==1);
print("# n<=N with e(BD)=2   : ", #select(t->t==2, vD), "    with e(BD)<=1 : ", #select(t->t<=1, vD));
print("# n<=N with e(Bnew)=2 : ", #select(t->t==2, vN), "    with e(Bnew)<=1 : ", #select(t->t<=1, vN));
print("");
print("e(BD_n)   n=1..40 : ", vector(40,i,vD[i]));
print("e(Bnew_n) n=1..40 : ", vector(40,i,vN[i]));
print("");
LD = select(t->vD[t]==2, vector(NA,i,i)); LN = select(t->vN[t]==2, vector(NA,i,i));
print("largest 15 n<=N with e(BD)=2   : ", vector(15,i,LD[#LD-15+i]));
print("largest 15 n<=N with e(Bnew)=2 : ", vector(15,i,LN[#LN-15+i]));
print("");
print("longest gap between consecutive n with e(BD)=2   : ", vecmax(vector(#LD-1,i,LD[i+1]-LD[i])));
print("longest gap between consecutive n with e(Bnew)=2 : ", vecmax(vector(#LN-1,i,LN[i+1]-LN[i])));
print("all A_n in Z ? ", vecmax(vector(NA+1,i,denominator(A[i])))==1);
quit;
