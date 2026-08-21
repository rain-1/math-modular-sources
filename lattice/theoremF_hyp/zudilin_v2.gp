/* zudilin_v2.gp -- exact 2-adic valuation law for Zudilin's Catalan row.
   Claim (verified n<=1200):  v_2(a_n) = -4n + 2*s_2(n).                   */
read("lattice/euler_criterion/rows.gp");
M=1200; Z=rowZud(M)[1]; ok=1;
for(n=1,M, if(valuation(Z[n+1],2) != -4*n+2*vecsum(digits(n,2)), ok=0; print("FAIL n=",n)));
print("v_2(a_n) = -4n + 2 s_2(n) for 1<=n<=", M, " : ", if(ok,"TRUE","FALSE"));
/* same for b_n */
P=rowZud(M)[2]; ok2=1; bad=List();
for(n=1,M, my(d=valuation(P[n+1],2)+4*n-2*vecsum(digits(n,2))); if(d!=0 && #bad<8, listput(bad,[n,d])));
print("v_2(b_n) - (-4n+2 s_2(n)) first deviations: ", Vec(bad));
quit;
