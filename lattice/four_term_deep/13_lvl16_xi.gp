default(parisize, 8000000000);
default(realprecision, 120);
L  = readvec("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_A.txt");
LB = readvec("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_B.txt");
NA = #L-1;
/* recurrence coefficient polynomials in nn, index s=0..r for u_{n+1-s} */
{Q6 = [ nn^2+2*nn+1, 14*nn^2+2*nn, 80*nn^2-120*nn+56, 240*nn^2-720*nn+576, 384*nn^2-1728*nn+1984, 256*nn^2-1536*nn+2304 ];}
{Q5 = [ nn^3+2*nn^2+nn, 10*nn^3+6*nn^2, 40*nn^3-24*nn^2-16*nn+16, 80*nn^3-144*nn^2-32*nn+128, 64*nn^3-192*nn^2+256 ];}
{run(QQ, u1, NMX) =
  my(r=#QQ-1, U=vector(NMX+1));
  U[1] = if(u1==0, 1, 0); U[2] = u1;
  for(n=1, NMX-1,
     my(s = 0);
     for(j=1, r, my(idx = n+1-j); if(idx>=0, s += subst(QQ[j+1],nn,n)*U[idx+1]));
     U[n+2] = -s/subst(QQ[1],nn,n));
  U;}
NMX = 1100;
A6 = run(Q6, 0, NMX); W6 = run(Q6, 1, NMX);
A5 = run(Q5, 0, NMX); W5 = run(Q5, 1, NMX);
print("a_n from 6-term rec matches file (n<=",NA,")? ", vector(NA+1,i,A6[i])==L);
print("a_n from 5-term cubic rec matches file? ", vector(NA+1,i,A5[i])==L);
print("w_n (6-term) first 8: ", vector(8,i,W6[i]));
print("w_n (5-term cubic) first 8: ", vector(8,i,W5[i]));
print("W6 == W5 up to n=",NMX,"? ", W6==W5);
{ for(k=1,6, my(m=NMX-5*(k-1));
   print("  xi6(",m,") = ", W6[m+1]/A6[m+1]*1.0)); }
xi6 = W6[NMX+1]/A6[NMX+1]*1.0;
print("diagnostic |xi6(",NMX,")-xi6(",NMX-1,")| = ", abs(xi6 - W6[NMX]/A6[NMX]*1.0));
xi5 = W5[NMX+1]/A5[NMX+1]*1.0;
print("xi5 = ", xi5, "   |xi5-xi5prev| = ", abs(xi5 - W5[NMX]/A5[NMX]*1.0));
print("-G/2 = ", -Catalan/2);
print("xi6 - (-G/2) = ", xi6 + Catalan/2);
print("xi6 / G = ", xi6/Catalan);
dn(n) = lcm(vector(max(n,1),j,j));
{ for(kx=0,6, my(ok=1,bad=0); for(m=1,60, if(denominator(dn(m)^kx*W6[m+1])!=1, ok=0; bad=m; break));
   print("  W6: d_n^",kx," w_n integral n<=60? ", ok, if(bad,Str(" first fail n=",bad),""))); }
{ for(kx=0,6, my(ok=1,bad=0); for(m=1,60, if(denominator(dn(m)^kx*W5[m+1])!=1, ok=0; bad=m; break));
   print("  W5: d_n^",kx," w_n integral n<=60? ", ok, if(bad,Str(" first fail n=",bad),""))); }
print("modular b_n/a_n (n=",NA,") = ", LB[NA+1]/L[NA+1]*1.0);
FN = "/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_limits.txt";
write1(FN,"");
write(FN, Str("xi_canon_6term ", xi6));
write(FN, Str("xi_canon_5term_cubic ", xi5));
/* D-form check */
Rc = 1+14*t+80*t^2+240*t^3+384*t^4+256*t^5;
Sc = 1+16*t+120*t^2+480*t^3+960*t^4+768*t^5;
Vc = 16*t^2+96*t^3+256*t^4+256*t^5;
{ my(ok=1); for(j=0,5, my(pj = polcoeff(Rc,j)*(nn+1-j)*(nn-j)+polcoeff(Sc,j)*(nn+1-j)+polcoeff(Vc,j));
    if(pj != Q6[j+1], ok=0; print("  MISMATCH j=",j,": ",pj," vs ",Q6[j+1])));
  print("D-form (Rc,Sc,Vc) reproduces the six-term recurrence? ", ok); }
Tc = Sc - t*deriv(Rc,t);
print("Rc factor: ", factor(Rc));  print("Sc factor: ", factor(Sc));
print("Vc factor: ", factor(Vc));  print("Tc factor: ", factor(Tc));
/* exponents at the double root t=-1/4, exactly */
{ my(u='u, t0=-1/4, p2=subst(t^2*Rc,t,t0+u), p1=subst(t*Sc,t,t0+u), p0=subst(Vc,t,t0+u));
  print("at t=-1/4: val_u(t^2 Rc)=",valuation(p2,u)," val_u(t Sc)=",valuation(p1,u)," val_u(Vc)=",valuation(p0,u));
  my(AA=polcoeff(p2,2,u), BB=polcoeff(p1,1,u), CC=polcoeff(p0,0,u));
  print("   A=",AA," B=",BB," C=",CC,"   indicial: ",AA*rr*(rr-1)+BB*rr+CC,"  roots ",polroots(AA*rr*(rr-1)+BB*rr+CC)); }
/* rho at the three simple roots, exactly via resultants */
{ foreach([2*t+1, 8*t^2+4*t+1], fac,
    print("  factor ",fac,": rho = -Tc/(t Rc') mod fac -> ",
       lift(Mod(-Tc,fac)*Mod(t*deriv(Rc,t),fac)^(-1))) ); }
{ my(r5=polcoeff(Rc,5), s5=polcoeff(Sc,5), v5=polcoeff(Vc,5));
  print("t=inf indicial (y~t^-mu): ", r5*mu*(mu+1)-s5*mu+v5, "  roots ", polroots(r5*mu*(mu+1)-s5*mu+v5)); }
print("P_5(n) = ", Q6[6], " = ", factor(Q6[6]));
print("P_4(n) cubic row = ", Q5[5], " = ", factor(Q5[5]));
print("primitivity: mu^s | P_s for mu=2? ", vector(5,s,denominator(polcoeff(Q6[s+1],0)/2^s)==1 && denominator(polcoeff(Q6[s+1],1)/2^s)==1 && denominator(polcoeff(Q6[s+1],2)/2^s)==1));
{ for(mu=2,4, my(ok=1); for(nx=0,NA, if(denominator(L[nx+1]/mu^nx)!=1, ok=0; break));
    print("  a_n/",mu,"^n integral all n<=",NA,"? ", ok)); }
{ my(bad=0); for(nx=0,NA, if(denominator(L[nx+1])!=1, bad=nx;break)); print("  a_n integral to n=",NA,"? ",bad==0); }
quit;
