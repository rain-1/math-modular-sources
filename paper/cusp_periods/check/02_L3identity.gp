/* 02_L3identity.gp -- Re L(3,psi_4) = phi^5 Im L(3,psi_4), 60 digits, and the exact
   generalised-Bernoulli derivation (Proposition 5.3 of the paper).            */
default(realprecision, 70);
G5 = znstar(5,1); Lf = lfuncreate([G5,[1]]);           /* psi_4 with psi_4(2) = i */
print("psi_4(2) = ", chareval(G5,[1],2,[I,4]));
L3 = lfun(Lf, 3); ph = (1+sqrt(5))/2; ph5 = ph^5;
print("L(3,psi_4) = ", L3);
print("Re/Im - phi^5 = ", real(L3)/imag(L3) - ph5);
/* exact formula: L(3,chi) = -i (2 pi)^3 tau(chi) B_{3,chibar} / (12 N^3) */
B3(x) = x^3 - 3/2*x^2 + x/2;
ps = [1,I,-I,-1,0]; psb = [1,-I,I,-1,0];
Bg(v) = 25*sum(a=1,5, v[a]*B3(a/5));
tau(v) = sum(a=1,5, v[a]*exp(2*Pi*I*a/5));
Lpred = -I*(2*Pi)^3*tau(ps)*Bg(psb)/(12*125);
print("B_{3,psibar_4} = ", Bg(psb), "   tau(psi_4) = ", tau(ps));
print("Bernoulli formula - lfun = ", Lpred - L3);
s1 = sin(2*Pi/5); s2 = sin(Pi/5);
print("Re L3 - (16 pi^3/7500)(12 s1 + 6 s2) = ", real(L3) - 16*Pi^3/7500*(12*s1+6*s2));
print("Im L3 - (16 pi^3/7500)(12 s2 - 6 s1) = ", imag(L3) - 16*Pi^3/7500*(12*s2-6*s1));
print("s1/s2 - phi = ", s1/s2 - ph, "   (2 phi+1)/(2-phi) - phi^5 = ", (2*ph+1)/(2-ph) - ph5);
print("polar coefficient of Phi_new at cusp 0:   ", (1+I*ph5)*L3 + (1-I*ph5)*conj(L3));
print("polar coefficient of Phi_new at cusp 1/2: ", ((1+I*ph5)*I*L3 + (1-I*ph5)*(-I)*conj(L3))/8);
L2 = lfun(Lf, 2);
print("L(2,psi_4) = ", L2);
print("xi        = ", ph5*imag(L2) - real(L2));
print("sigma(xi) = ", -real(L2) - imag(L2)/ph5);
quit;
