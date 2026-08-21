/* ===========================================================================
   zeta5_check5.gp -- hypothesis (b) of Theorem F (DESCENT) for the zeta(5)
   hosts, and the 2-adic location of the branch points (folds).

   For a host coordinate t that is invariant under an Atkin-Lehner involution
   W_Q, the row descends to P^1_t only if the source Phi is a W_Q-eigenform
   with eigenvalue +1 (then F=Phi/(Dt) is W_Q-invariant of weight w, and
   H_{xi*}=F(Theta-xi*) is Gamma_t-invariant).

   Exact test.  Phi = sum_d c_d E_{w+2}(d tau) on level N.  Then
       (Phi |_{w+2} W_N)  =  sum_{d'} c'_{d'} E_{w+2}(d' tau),
       c'_{d'} = c_{N/d'} * d'^{w+2} / N^{(w+2)/2}.
   =========================================================================== */
\p 40

frick(N, dv, cv, k) =
{ my(r=vector(#dv), j);
  for(i=1,#dv, my(dp=dv[i], d=N/dp); j=0;
     for(m=1,#dv, if(dv[m]==d, j=m));
     if(j==0, error("divisor set not W_N-stable"));
     r[i] = cv[j]*dp^k/N^(k/2));
  r;
}

test(name, N, dv, cv) =
{ my(w = frick(N,dv,cv,6), ev);
  print("--- ", name, "  (N=",N,")");
  print("    c        = ", cv);
  print("    c|W_",N,"  = ", w);
  ev = w[1]/cv[1];
  if(w == ev*cv,
     print("    ==> Phi IS a W_",N," eigenform, eigenvalue ", ev,
           "   -> descent onto a W_",N,"-invariant t: ", if(ev==1,"OK","FAILS (eigenvalue -1)")),
     print("    ==> Phi is NOT a W_",N," eigenform  -> HYPOTHESIS (b) FAILS"));
}

print("############ hypothesis (b): is the source a Fricke eigenform? ############");
d12=[1,2,3,4,6,12]; d16=[1,2,4,8,16];
test("zeta(5) level-16 source  L(Phi,5)=(217/1024)zeta(5)", 16, d16, [1,-85,1428,-5440,4096]);
test("zeta(5) level-12 source  L(Phi,5)=(31/192)zeta(5)  [the puzzle]", 12, d12, [1,-113,567,112,-1863,1296]);
test("zeta(5) level-12 archive  (25/144), Domb", 12, d12, [1,-104,351,832,-2808,1728]);
test("zeta(5) level-12 archive  (11/144), anti-Fricke", 12, d12, [1,-176,2079,-4928,4752,-1728]);
print("");

print("############ numerical confirmation in the upper half plane ############");
etaq(L,tau) = { my(r=1.0); for(i=1,#L, r *= eta(L[i][1]*tau,1)^L[i][2]); r; }
E6s(tau) = { my(NT=400, qq=exp(2*Pi*I*tau)); 1 - 504*sum(n=1,NT, sigma(n,5)*qq^n); }
PHIv(cv,dv,tau) = sum(i=1,#dv, cv[i]*(1-E6s(dv[i]*tau))/504);
tau = 0.31 + 0.77*I;

Lh = [[1,3],[2,-2],[3,-1],[4,1],[6,2],[12,-3]];
h  = etaq(Lh,tau); hW = etaq(Lh,-1/(12*tau));
print("level-12 h_12 host:  h*(h|W_12) = ", h*hW, "   (=12  =>  t=h/((h+3)(h+4)) is W_12-invariant)");
print("   t(W_12 tau)/t(tau) = ", (hW/((hW+3)*(hW+4)))/(h/((h+3)*(h+4))));
print("   (Phi_12|_6 W_12)/Phi_12 = ", (12^-3*tau^-6*PHIv([1,-113,567,112,-1863,1296],d12,-1/(12*tau)))/PHIv([1,-113,567,112,-1863,1296],d12,tau));
print("      -> not +-1: H_{xi*} is NOT Gamma_t-invariant.  (b) FAILS.");
Lx = [[1,-2],[2,1],[8,-1],[16,2]];
x = etaq(Lx,tau); xW = etaq(Lx,-1/(16*tau));
print("level-16 host:  8*x*(x|W_16) = ", 8*x*xW, "   (=1  =>  t=x/(8x^2+2x+1) is W_16-invariant)");
print("   t(W_16 tau)/t(tau) = ", (xW/(8*xW^2+2*xW+1))/(x/(8*x^2+2*x+1)));
print("   (Phi_16|_6 W_16)/Phi_16 = ", (16^-3*tau^-6*PHIv([1,-85,1428,-5440,4096],d16,-1/(16*tau)))/PHIv([1,-85,1428,-5440,4096],d16,tau));
print("      -> +1 exactly: (b) HOLDS.");
print("");

print("############ where the folds (branch points) sit 2-adically ############");
print("level-12 h_12:  t = h/((h+3)(h+4)) ; dt/dh = 0 at h^2 = 12, i.e. h = +-2sqrt3");
print("   t(2sqrt3) = 7 - 4*sqrt(3) = ", 7-4*sqrt(3), " ;  t(-2sqrt3) = 7 + 4*sqrt(3)");
print("   minimal polynomial of 7-4sqrt3 : x^2 - 14x + 1 ;  over Q_2:");
print("     ", factorpadic('x^2-14*'x+1, 2, 20));
print("     Newton polygon (0,0)-(2,0): BOTH roots are 2-adic UNITS, v_2 = 0, |t|_2 = 1.");
print("   ==> the fold sits ON the unit circle |t|_2 = 1: rho = 1, sigma_2 = 0.");
print("");
print("level-16:  t = w/(1+2w), w = x/(1+8x^2) ; dw/dx = 0 at x = +-1/(2sqrt2)");
print("   folds t = (-1 +- 2sqrt2)/14 ; minimal polynomial 28t^2 + 4t - 1 ; over Q_2:");
print("     ", factorpadic(28*'x^2+4*'x-1, 2, 20));
print("   v_2 of its roots = ", (valuation(-1,2)-valuation(28,2))/2, "   (Newton polygon (0,0)-(2,2) slope -1)");
print("   other singular points t = -1/4, -1/2 : v_2 = -2, -1.");
print("   ==> the 2-adically NEAREST singular t has v_2 = -1 (cusp -1/2 and the folds),");
print("       so rho = min|s|_2 = 2^1 = 2 and the PREDICTED SLOPE is sigma_2 = log_2 rho = 1.");
print("");
print("############ WHY (a) and (b) are incompatible at level 12 ############");
crit2(c) = [c[1]+c[2]+c[4], c[3]+c[5]+c[6]];
cp=[1,-104,351,832,-2808,1728]; cm=[1,-176,2079,-4928,4752,-1728]; cx=[1,-113,567,112,-1863,1296];
print("  Fricke-EVEN (25/144): criterion sums = ", crit2(cp), "   passes (a)? ", crit2(cp)==[0,0]);
print("  Fricke-ODD  (11/144): criterion sums = ", crit2(cm), "   passes (a)? ", crit2(cm)==[0,0]);
print("  the puzzle  (31/192): criterion sums = ", crit2(cx), "   passes (a)? ", crit2(cx)==[0,0]);
mix = (7*cp + cm)/8;
print("  (7*Phi_even + Phi_odd)/8 = ", mix, "   equals the 31/192 source? ", mix==cx);
print("  L-values: (7/8)(25/144) + (1/8)(11/144) = ", 7/8*25/144 + 1/8*11/144);
print("  ==> on the 2-dimensional level-12 purified span {Phi_even, Phi_odd}, the Euler-factor");
print("      criterion (a) is a codimension-2 linear condition whose unique solution line is a");
print("      NON-EIGEN mixture.  (a) and (b) are mutually exclusive at N=12: any source passing");
print("      (a) fails to descend to a Fricke-invariant Hauptmodul.");
print("DONE"); quit;
