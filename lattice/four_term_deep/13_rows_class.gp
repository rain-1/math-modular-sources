/* Uniform treatment: class invariants + canonical companion + xi, for the three rows. */
default(parisize, 8000000000);
default(realprecision, 150);
OUT="/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/";
dn(n) = lcm(vector(max(n,1),j,j));
{L16 = [nn^2+2*nn+1, 14*nn^2+2*nn, 80*nn^2-120*nn+56, 240*nn^2-720*nn+576, 384*nn^2-1728*nn+1984, 256*nn^2-1536*nn+2304];}
{P3 = [nn^2+2*nn+1, -5*nn^2-5*nn-2, -4*nn^2+16*nn-4, 40*nn^2-120*nn+76, -16*nn^2+32*nn+32, -80*nn^2+560*nn-992, 64*nn^2-512*nn+1024];}
{P6 = [nn^2+2*nn+1, 9*nn^2+nn, 20*nn^2-40*nn+20, -20*nn^2+20*nn+4, -96*nn^2+352*nn-336, -64*nn^2+320*nn-400];}
{runrec(QQ, u1, NMX) =
  my(r=#QQ-1, U=vector(NMX+1));
  U[1] = if(u1==0,1,0); U[2] = u1;
  for(n=1,NMX-1, my(s=0);
     for(j=1,r, my(idx=n+1-j); if(idx>=0, s += subst(QQ[j+1],nn,n)*U[idx+1]));
     U[n+2] = -s/subst(QQ[1],nn,n));
  U;}
{classof(QQ, tag) =
  my(r=#QQ-1);
  my(Rc=0, Sc=0, Vc=0);
  for(j=0,r,
    my(rj = polcoeff(QQ[j+1],2,nn), sj = polcoeff(QQ[j+1],1,nn) - rj*(1-2*j));
    my(vj = polcoeff(QQ[j+1],0,nn) - rj*(j^2-j) - sj*(1-j));
    Rc += rj*t^j; Sc += sj*t^j; Vc += vj*t^j);
  print("\n##### ", tag);
  my(ok=1); for(j=0,r, if(polcoeff(Rc,j)*(nn+1-j)*(nn-j)+polcoeff(Sc,j)*(nn+1-j)+polcoeff(Vc,j) != QQ[j+1], ok=0));
  print("  D-form check: ", ok);
  print("  Rc = ", Rc, "  = ", factor(Rc));
  print("  Sc = ", Sc);
  print("  Vc = ", Vc, if(Vc==0,"",Str("  = ",factor(Vc))));
  my(Tc = Sc - t*deriv(Rc,t));
  print("  Tc = ", Tc);
  print("  exponents at t=0: roots of ", polcoeff(Rc,0),"*rr*(rr-1)+",polcoeff(Sc,0),"*rr+",polcoeff(Vc,0), " -> ", polroots(polcoeff(Rc,0)*rr*(rr-1)+polcoeff(Sc,0)*rr+polcoeff(Vc,0)));
  my(fa = factor(Rc));
  for(i=1,matsize(fa)[1], my(fac=fa[i,1], mult=fa[i,2]);
    if(mult==1,
      print("  factor ", fac, " (simple):  rho = -Tc/(t Rc') = ", lift(Mod(-Tc,fac)*Mod(t*deriv(Rc,t),fac)^(-1)), "   [exponents {0, rho}]"),
      /* multiple root: indicial via local expansion at each root */
      my(rts = polroots(fac));
      for(z=1,#rts, my(t0=rts[z], u='u);
        my(p2 = subst(t^2*Rc,t,t0+u), p1 = subst(t*Sc,t,t0+u), p0 = subst(Vc,t,t0+u));
        my(vv = 10^10);
        for(e=0,10, if(abs(polcoeff(p2,e,u))>1e-60, vv=min(vv,e-2); break));
        for(e=0,10, if(abs(polcoeff(p1,e,u))>1e-60, vv=min(vv,e-1); break));
        for(e=0,10, if(abs(polcoeff(p0,e,u))>1e-60, vv=min(vv,e); break));
        my(AA=polcoeff(p2,vv+2,u), BB=polcoeff(p1,vv+1,u), CC=polcoeff(p0,vv,u));
        print("  root ", t0, " (mult ",mult,"): indicial ",AA,"*rr(rr-1)+",BB,"*rr+",CC,"  -> exponents ", polroots(AA*rr*(rr-1)+BB*rr+CC)))));
  my(muP = subst(QQ[r+1], nn, r-1-mu));
  print("  P_",r,"(n) = ", QQ[r+1], " = ", factor(QQ[r+1]));
  print("  exponents at infinity (y ~ t^-mu): roots of P_",r,"(",r-1,"-mu) = ", muP, " -> ", polroots(muP));
  my(mr = polroots(muP)); print("  delta_infinity = |mu_2-mu_1| = ", abs(mr[2]-mr[1]));
  my(mods = vecsort(vector(poldegree(Rc), i, abs(polroots(Rc)[i]))));
  print("  |t| sorted: ", vector(#mods,i,mods[i]));
  print("  lambda_1 = ", 1/mods[1], "   lambda_2 (next DISTINCT modulus) = ",
     my(d=0); for(i=2,#mods, if(abs(mods[i]-mods[1])>1e-40, d=1/mods[i]; break)); d);
  [Rc,Sc,Vc];}
{doit(QQ, afile, tag, NMX) =
  classof(QQ, tag);
  my(A = readvec(afile), NA=#A-1);
  my(AA = runrec(QQ,0,NMX), W = runrec(QQ,1,NMX));
  print("  recurrence reproduces the modular a_n (n<=",min(NA,NMX),")? ", vector(min(NA,NMX)+1,i,AA[i])==vector(min(NA,NMX)+1,i,A[i]));
  print("  w_n first 8: ", vector(8,i,W[i]));
  my(xi = W[NMX+1]/AA[NMX+1]*1.0, xip = W[NMX]/AA[NMX]*1.0);
  print("  xi_canonical(n=",NMX,") = ", xi);
  print("  |xi(n)-xi(n-1)| = ", abs(xi-xip));
  for(kx=0,5, my(ok=1); for(m=1,60, if(denominator(dn(m)^kx*W[m+1])!=1, ok=0;break));
     if(ok, print("  k_den(w) = ", kx); break));
  my(bad=0); for(nx=0,NA, if(denominator(A[nx+1])!=1, bad=nx;break));
  print("  a_n integral to n=",NA,"? ", bad==0);
  for(mu=2,4, my(o=1); for(nx=0,NA, if(denominator(A[nx+1]/mu^nx)!=1, o=0;break)); if(o, print("  *** a_n/",mu,"^n integral -- NOT primitive")));
  print("  content gcd(a_0..a_",NA,") = ", gcd(vector(NA+1,i,A[i])));
  xi;}
x16 = doit(L16, Str(OUT,"lvl16_A.txt"), "level 16 (six-term)", 900);
x3  = doit(P3,  Str(OUT,"lvl12_p3_A.txt"), "Gamma_0(12) placement 3, c=3 (seven-term)", 500);
x6  = doit(P6,  Str(OUT,"lvl12_p6_A.txt"), "Gamma_0(12) placement 6, c=6 (six-term)", 500);
FN = Str(OUT,"lvl12_limits.txt"); write1(FN,"");
write(FN, Str("xi_canon_lvl12_p3 ", x3));
write(FN, Str("xi_canon_lvl12_p6 ", x6));
print("\ncontrols: -G/27 = ", -Catalan/27, "   G/9 = ", Catalan/9);
print("          -(4/9)zeta(2) = ", -4*zeta(2)/9, "   (2/9)zeta(2) = ", 2*zeta(2)/9);
print("x3 vs -G/27 : ", x3 + Catalan/27, "    x3 vs -(4/9)z2: ", x3 + 4*zeta(2)/9);
print("x6 vs G/9   : ", x6 - Catalan/9,  "    x6 vs (2/9)z2 : ", x6 - 2*zeta(2)/9);
quit;
