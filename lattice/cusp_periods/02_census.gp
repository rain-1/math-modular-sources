/* 02_census.gp -- the cusp-period census: every cusp of Gamma_0(N_Phi) (resp.
   Gamma_1(5)) for each of the twelve table sources and the two new Gamma_1(5)
   directions.  Prints Pi, the polar coefficient rho, the fold-regularity flag,
   and the orientation factor Pi/Pi(cusp 0).                                  */
read("lib.gp");
default(realprecision, 60);
TOL = 10^(-40);
fmt(z) = if(abs(imag(z)) < TOL, Str(real(z)), Str(real(z)," + ",imag(z),"*I"));
ratapp(z) = my(b=bestappr(z, 10^20)); if(abs(b-z) < TOL, Str(b), "-");
docusp(s, a, c, P0) = my(b=cperiod(s,a,c)); my(fr=if(abs(b[2])<TOL,"YES","no ")); my(of=if(abs(P0)>TOL, ratapp(b[1]/P0), "-")); printf("  %2d/%-3d  foldreg %s  Pi = %s\n", a, c, fr, fmt(b[1])); printf("            rho = %s   Pi/Pi(0) = %s\n", fmt(b[2]), of);
dorow(s, cl) = my(P0=cperiod(s,0,1)[1]); print("== ", s[1], "   r = ", s[2], "   level ", s[3]); print("   L(Phi,r) = Pi(cusp 0) = ", fmt(P0)); for(i=1,#cl, docusp(s, cl[i][1], cl[i][2], P0)); print("   Pi(cusp infinity) = 0 by construction"); print();
for(i=1,#SRC, dorow(SRC[i], cusplist(SRC[i][3])));
print("=== Gamma_1(5): the two quartic directions (complex combinations)");
ph5 = (11+5*sqrt(5))/2;
Ep = srcbyname("E3ps"); Em = srcbyname("E3psb");
mixp(cf1, cf2, a, c) = my(u=cperiod(Ep,a,c)); my(v=cperiod(Em,a,c)); [cf1*u[1]+cf2*v[1], cf1*u[2]+cf2*v[2]];
G15 = [[0,1],[1,2],[2,5],[1,5]];
shownew(nm, c1, c2) = print("-- ", nm); for(i=1,#G15, my(cu=G15[i]); my(b=mixp(c1,c2,cu[1],cu[2])); printf("   %d/%d  foldreg %s   Pi = %s   rho = %s\n", cu[1], cu[2], if(abs(b[2])<TOL,"YES","no "), fmt(b[1]), fmt(b[2])));
shownew("Phi_new  = (1+i phi^5) E^{psi4,1} + (1-i phi^5) E^{psibar4,1}", 1+I*ph5, 1-I*ph5);
shownew("Phi'_new = (1-i/phi^5) E^{psi4,1} + (1+i/phi^5) E^{psibar4,1}", 1-I/ph5, 1+I/ph5);
shownew("R3 = E^{psi4,1}+E^{psibar4,1}", 1, 1);
shownew("R4 = i(E^{psi4,1}-E^{psibar4,1})", I, -I);
quit;
