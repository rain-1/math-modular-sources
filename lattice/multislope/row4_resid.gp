default(parisizemax, 6000000000);
default(realprecision,120);
An = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An.txt");
Bn = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn.txt");
ciB = read("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_recB.txt");
QB = vector(36, k, my(i=k-1); sum(j=0,7, ciB[i*8+j+1]*'n^j));
resid(x,n) = { my(sm=0); for(i=0,35, if(n-i>=0, sm += subst(QB[i+1],'n,n)*x[n-i+1])); sm }
cn = vector(121, k, resid(Bn, k-1));
print("c_n = L[B]_n for n=0..10: ", vector(11,k,cn[k]));
print("ratios c_n/c_{n-1} for n=2..8: ", vector(7,k, my(n=k+1); cn[n+1]/cn[n]));
print("factored ratios: ", vector(7,k, my(n=k+1); [n, factor(cn[n+1]/cn[n])~[1,]]));
\\ fit first-order recurrence u(n) c_n + v(n) c_{n-1} = 0 with deg <= DD
pp = 2^61-1;
cm = vector(121, k, if(cn[k]==0, Mod(0,pp), Mod(numerator(cn[k]),pp)/Mod(denominator(cn[k]),pp)));
kd(DD) = {
  my(nc = 2*(DD+1));
  my(mat = matrix(115-DD, nc, a, b, cm[(DD+5)+a-1-((b-1)\(DD+1))+1] * Mod((DD+5)+a-1,pp)^((b-1)%(DD+1))));
  nc - matrank(mat)
}
print("first-order fit to c_n: kernel dims for deg D=0..12: ", vector(13,j,kd(j-1)));
print("");
print("--- archimedean behaviour of the companions ---");
ciBv = ciB;
NC=260;
companion(s) = { my(x=vector(NC+1)); x[s+1]=1; for(n=s+1,NC, my(sm=0); for(i=1,35, if(n-i>=0, sm += subst(QB[i+1],'n,n)*x[n-i+1])); x[n+1] = -sm/n^7); x }
X1=companion(1); X2=companion(2);
print("X1_n/A_n at n=200,230,260: ", vector(3,k,my(n=170+30*k); 1.0*X1[n+1]/An[n+1]));
print("X2_n/A_n at n=200,230,260: ", vector(3,k,my(n=170+30*k); 1.0*X2[n+1]/An[n+1]));
print("B_n/A_n  at n=260: ", 1.0*Bn[261]/An[261]);
print("successive diffs |X1_n/A_n - X1_{n-1}/A_{n-1}|^(1/n) at n=260: ", abs(1.0*(X1[261]/An[261]-X1[260]/An[260]))^(1.0/260));
quit
