default(parisizemax, 6000000000);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
ciB = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt");
ciA = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recA.txt");
mk(ci,r,D) = vector(r+1, k, my(i=k-1); sum(j=0,D, ci[i*(D+1)+j+1]*'n^j));
QB = mk(ciB,35,7); QA = mk(ciA,22,23);
resid(Q,r,x,n) = { my(sm=0); for(i=0,r, if(n-i>=0, sm += subst(Q[i+1],'n,n)*x[n-i+1])); sm }
print("residual L_B[B_n] for n=0..45:");
{my(v=vector(46,k,resid(QB,35,Bn,k-1))); for(k=1,46, if(v[k]!=0, print("  n=",k-1," resid nonzero"), ))}
print("  count nonzero n=0..400: ", sum(k=0,400, if(resid(QB,35,Bn,k)!=0,1,0)));
print("  first 6 nonzero n: ", select(k->resid(QB,35,Bn,k)!=0, vector(401,k,k-1))[1..min(6,#select(k->resid(QB,35,Bn,k)!=0, vector(401,k,k-1)))]);
print("residual L_A[B_n], count nonzero n=0..400: ", sum(k=0,400, if(resid(QA,22,Bn,k)!=0,1,0)));
print("  nonzero n (first 8, shape A): ", select(k->resid(QA,22,Bn,k)!=0, vector(401,k,k-1))[1..8]);
print("residual L_A[A_n], count nonzero n=0..400: ", sum(k=0,400, if(resid(QA,22,An,k)!=0,1,0)));
quit
