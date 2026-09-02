\\ The complete HOLOMORPHIC fold-regular source space on each Fricke host level:
\\ V(N) = { Phi in M_4(Gamma_0(N)) : Phi|W_N = -Phi, a_0(Phi) = 0 },
\\ split as V = (Eisenstein part) + (cuspidal part).
default(parisizemax, 16000000000);
read("lib.gp");
doN(N) = my(mf,d,Bb,ai,M,Km,V,rows,W,g,cv,mfc,dc,ns); mf=mfinit([N,4],4); d=mfdim(mf); Bb=mfbasis(mf); ai=mfatkininit(mf,N); M=matrix(d,d); for(j=1,d, my(vv=mftobasis(mf,mfatkin(ai,Bb[j]))); for(t=1,d, M[t,j]=vv[t])); Km=matker(M+1); mfc=mfinit([N,4],1); dc=mfdim(mfc); if(matsize(Km)[2]==0, print("N=",N,"  dimM4=",d," dimS4=",dc,"  dim(W=-1)=0   dim V = 0"); return(0)); rows=vector(matsize(Km)[2], j, mfcoefs(mflinear(mf,Km[,j]),0)[1]); W=matker(Mat(rows)); V=Km*W; ns=0; if(dc>0, my(Bc=mfbasis(mfc), aic=mfatkininit(mfc,N), Mc=matrix(dc,dc)); for(j=1,dc, my(vv=mftobasis(mfc,mfatkin(aic,Bc[j]))); for(t=1,dc, Mc[t,j]=vv[t])); ns=matsize(matker(Mc+1))[2]); print("N=", N, "  dimM4=", d, "  dimS4=", dc, "   dim(W_N=-1 on M4) = ", matsize(Km)[2], "   dim V = ", matsize(V)[2], "   [ V_Eis = ", matsize(V)[2]-ns, " , V_cusp = ", ns, " ]"); for(j=1,matsize(V)[2], g=mflinear(mf, primvec(V[,j]~)~); print("   V[",j,"] = ", mfcoefs(g,14))); matsize(V)[2];
lev = [5,6,7,8,9,10,12,18];
for(i=1,#lev, doN(lev[i]); print(""));
quit;
