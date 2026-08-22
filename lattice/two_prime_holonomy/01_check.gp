default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/01_lib.gp");
PX = 40; PY = 18;
HA = solveode(0,0,PX); HB = solveode(1,0,PX); HC = solveode(0,1,PX);
print("HA[0..3] = ", vector(4,i,HA[i]));
print("HB[0..3] = ", vector(4,i,HB[i]));
print("HC[0..3] = ", vector(4,i,HC[i]));
\\ Li2
D2 = lcm(vector(PX,i,i))^2;
Li2 = vector(PX); for(n=1,PX-1, Li2[n+1] = D2/n^2);
Li2w = composew(Li2, PX);
dif = vector(PX, i, Li2[i]-Li2w[i]);
pr4 = antimul(dif, PX);
B4i = toy(pr4, PX, PY);
B4 = vector(PY, i, B4i[i]/D2);
print("B4[1..3] = ", vector(3,i,B4[i+1]));
\\ w-series sanity
tst = vector(PX); tst[2]=1;
print("w[0..4] = ", vector(5,i,composew(tst,PX)[i]));
\\ G sanity for (1,-3,5)
H = vector(PX, i, HA[i]-3*HB[i]+5*HC[i]);
DD = lcm(vector(PX,i,denominator(H[i])));
Hi = vector(PX, i, H[i]*DD);
Gx = vector(PX, i, Hi[i]+composew(Hi,PX)[i]);
Gi = toy(Gx, PX, PY);
G = vector(PY, i, Gi[i]/DD);
print("G[0..4] = ", vector(5,i,G[i]));
