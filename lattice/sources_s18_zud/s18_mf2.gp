read("/home/ubuntu/code/math-modular-sources/lattice/sources_s18_zud/s18_setup.gp");
Z = x*deriv(U18/x)/(U18/x);
W = -F18*Z;
nc = 40;
wv = vector(nc+1, i, polcoeff(W, i-1));
mf4E = mfinit([18,4],3); mf4S = mfinit([18,4],1);
BE = mfbasis(mf4E); BS = mfbasis(mf4S);
print("dim Eis=",#BE,"  dim cusp=",#BS);
cols = concat(vector(#BE,i,mfcoefs(BE[i],nc)~), vector(#BS,i,mfcoefs(BS[i],nc)~));
Mx = Mat(cols);
sol = matsolve(Mx[1..#cols,], wv[1..#cols]~);
print("solution (Eis then cusp): ", sol~);
print("check residual: ", Mx*sol - wv~);
Epart = sum(i=1,#BE, sol[i]*mfcoefs(BE[i],nc));
Cpart = sum(i=1,#BS, sol[#BE+i]*mfcoefs(BS[i],nc));
print("Eisenstein part of W: ", Epart);
print("cuspidal part of W  : ", Cpart);
{ for(i=1,#BE, print("  Eis basis ",i,": ", mfcoefs(BE[i],12), "  params ", mfparams(BE[i]))); }
{ for(i=1,#BS, print("  Cusp basis ",i,": ", mfcoefs(BS[i],12))); }
print("--- newform decomposition of S_4(18) ---");
S = mfinit([18,4],0);
print("new dim ", mfdim(S));
{ my(nf=mfeigenbasis(S)); for(i=1,#nf, print("  newform ",i,": ", mfcoefs(nf[i],14))); }
quit;
