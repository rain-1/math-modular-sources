/* 03_L3exact.gp -- the critical-L-value fold-regularity condition of (P1) at the
   cusp 1/2 of Gamma_1(5), and its EXACT algebraic proof.

   For an odd N-periodic h,  L(3,h) = i * sum_b hhat(b) * (2 pi)^3 B_3(b/N)/12,
   hhat(b) = (1/N) sum_a h(a) e(-ab/N), from  sum_{m>=1} sin(2 pi m x)/m^3 =
   (2 pi)^3 B_3(x)/12.  For a primitive odd chi mod N this is the generalized
   Bernoulli formula
        L(3,chi) = -i (2 pi)^3 tau(chi) B_{3,chibar} / (12 N^3).            */
default(realprecision, 80);
ee(x) = exp(2*Pi*I*x);
B3(x) = x^3 - 3/2*x^2 + 1/2*x;
ps4 = [1,I,-I,-1,0];
psb = [1,-I,I,-1,0];
pv5(v,n) = v[((n-1)%5)+1];
/* generalized Bernoulli B_{3,chi} = N^2 sum_a chi(a) B_3(a/N) */
Bgen(v) = 25*sum(a=1,5, pv5(v,a)*B3(a/5));
tau(v) = sum(a=1,5, pv5(v,a)*ee(a/5));
print("B_{3,psi4}     = ", Bgen(ps4));
print("B_{3,psibar4}  = ", Bgen(psb));
print("tau(psi4)      = ", tau(ps4));
Lpred = -I*(2*Pi)^3*tau(ps4)*Bgen(psb)/(12*125);
G5 = znstar(5,1);
Lact = lfun(lfuncreate([G5,[1]]),3);
print("L(3,psi4) predicted by Bernoulli/Gauss = ", Lpred);
print("L(3,psi4) from lfun                    = ", Lact);
print("difference                             = ", Lpred - Lact);
ph5 = (11+5*sqrt(5))/2;
print();
print("Re L(3,psi4) / Im L(3,psi4) = ", real(Lact)/imag(Lact));
print("phi^5                       = ", ph5);
print("difference                  = ", real(Lact)/imag(Lact) - ph5);
print();
print("EXACT: with s1 = sin(2 pi/5), s2 = sin(4 pi/5) = sin(pi/5):");
s1 = sin(2*Pi/5); s2 = sin(4*Pi/5);
print("  Re L(3,psi4) = (16 pi^3/7500)(12 s1 + 6 s2) = ", (16*Pi^3/7500)*(12*s1+6*s2));
print("  Im L(3,psi4) = (16 pi^3/7500)(12 s2 - 6 s1) = ", (16*Pi^3/7500)*(12*s2-6*s1));
print("  ratio = (2 s1 + s2)/(2 s2 - s1), and s1/s2 = 2 cos(pi/5) = phi :");
print("  s1/s2 = ", s1/s2, "   phi = ", (1+sqrt(5))/2);
print("  (2 phi + 1)/(2 - phi) = ", (2*(1+sqrt(5))/2+1)/(2-(1+sqrt(5))/2));
print("  phi^3 = 2 phi + 1 = ", ((1+sqrt(5))/2)^3, "   phi^-2 = 2 - phi = ", ((1+sqrt(5))/2)^(-2));
print("  hence the ratio is phi^3 * phi^2 = phi^5 EXACTLY.");
print();
print("=== the same statement as a vanishing of a generalized Bernoulli number");
print("v = 2 re(psi4) - 2 phi^5 im(psi4):  B_{3,v} = 2 B_{3,re} - 2 phi^5 B_{3,im}");
re4 = [1,0,0,-1,0]; im4 = [0,1,-1,0,0];
print("  B_{3,re4} = ", Bgen(re4), "   B_{3,im4} = ", Bgen(im4));
print("  ** the naive ratio B_{3,re4}/B_{3,im4} = ", Bgen(re4)/Bgen(im4), " is NOT phi^5:");
print("  the Gauss sum tau does not commute with re/im, so the L-ratio is not the B-ratio.");
print("  Correct exact statement: (1 + i phi^5) L(3,psi4) + (1 - i phi^5) L(3,psibar4) = 0,");
print("  value = ", (1+I*ph5)*Lact + (1-I*ph5)*conj(Lact));
print();
print("=== fold-regularity of Phi_new at cusp 1/2 vs cusp 0 (the (P1) L-condition)");
print("  Phi_new is inner (psi = quartic, phi = 1): polar coeff at a/c is");
print("  psi(c) c^{-r-1} L(r+1,psi) summed over the K-combination, r = 2:");
print("  at c = 1 (cusp 0):  L(3, (1+i phi^5) psi4 + (1-i phi^5) psibar4) = ", (1+I*ph5)*Lact+(1-I*ph5)*conj(Lact));
Lps = lfun(lfuncreate([G5,[1]]),3);
print("  at c = 2 (cusp 1/2): psi4(2) = i, psibar4(2) = -i, factor 2^{-3}:");
print("     (1/8)[ i(1+i phi^5) L(3,psi4) - i(1-i phi^5) L(3,psibar4) ] = ", (1/8)*(I*(1+I*ph5)*Lps - I*(1-I*ph5)*conj(Lps)));
quit;
