/* src.gp -- Fourier coefficients c(m) of the twelve Eisenstein sources,
   transcribed from paper/sections/02_sources.tex Table tab:sources
   (identical to lattice/thmB_exact/common.gp).                            */

sigin(chi,k,m) = if(m<1 || frac(m), 0, sumdiv(m, d, chi(m/d)*d^k));
sigout(chi,k,m) = if(m<1 || frac(m), 0, sumdiv(m, d, chi(d)*d^k));

ch1(n) = 1;
chm3(n) = kronecker(-3,n);
chm4(n) = kronecker(-4,n);
ch5(n) = kronecker(5,n);
nuD(d) = my(u=d%5); if(u==0, 0, [1,-2,2,-1][u]);

cmA(m) = sigout(chm3,2,m) - sigout(chm3,2,m/2);
cmB(m) = sigin(chm3,2,m) - 6*sigin(chm3,2,m/2) - 8*sigin(chm3,2,m/4);
cmC(m) = sigin(chm3,2,m) - 8*sigin(chm3,2,m/2);
cmD(m) = sumdiv(m, d, nuD(d)*d^2);
cmE(m) = sigin(chm4,2,m) - 8*sigin(chm4,2,m/2);
cmF(m) = sigin(chm3,2,m) - 7*sigin(chm3,2,m/2) - 8*sigin(chm3,2,m/4);
cmAl(m) = sigin(ch1,3,m) - 17*sigin(ch1,3,m/2) - 9*sigin(ch1,3,m/3) + 16*sigin(ch1,3,m/4) + 153*sigin(ch1,3,m/6) - 144*sigin(ch1,3,m/12);
cmGa(m) = sigin(ch1,3,m) - 28*sigin(ch1,3,m/2) + 63*sigin(ch1,3,m/3) - 36*sigin(ch1,3,m/6);
cmDe(m) = sigin(ch1,3,m) - 14*sigin(ch1,3,m/2) - sigin(ch1,3,m/3) + 16*sigin(ch1,3,m/4) + 14*sigin(ch1,3,m/6) - 16*sigin(ch1,3,m/12);
cmEp(m) = sigin(ch1,3,m) - 21*sigin(ch1,3,m/2) + 84*sigin(ch1,3,m/4) - 64*sigin(ch1,3,m/8);
cmZe(m) = sumdiv(m, d, chm3(d)*chm3(m/d)*d^3);
cmEt(m) = sigin(ch5,3,m) - 14*sigin(ch5,3,m/2) - 16*sigin(ch5,3,m/4);

{CMS = [cmA, cmB, cmC, cmD, cmE, cmF, cmAl, cmGa, cmDe, cmEp, cmZe, cmEt];}
