/* Fold-regularity + period + 2-adic profile for weight-3 chi_{-4} classes
   on the level-8 host (t,F8) and the level-16 host (x,F16).            */
read("00_setup.gp");
G = Catalan; Z2 = Pi^2/6;
H8  = mkhost(t8,F8);   A8v  = Aof(H8);
H16 = mkhost(x16,F16); A16v = Aof(H16);

/* class list: [name, base ("in"/"out"), P-list] */
{CL = [
 ["E_8   = (1-8V2)E            ", "in",  [[1,1],[2,-8]]],
 ["E_4   = (1-4V2)E            ", "in",  [[1,1],[2,-4]]],
 ["E_2   = (1-2V2)E            ", "in",  [[1,1],[2,-2]]],
 ["E     = E                   ", "in",  [[1,1]]],
 ["Psi   = (1-8V2)(1-4V2)E     ", "in",  [[1,1],[2,-12],[4,32]]],
 ["E_8'  = V2(1-8V2)E          ", "in",  [[2,1],[4,-8]]],
 ["A0    = (1-2V2)E   (notes)  ", "in",  [[1,1],[2,-2]]],
 ["T     = T                   ", "out", [[1,1]]],
 ["T4    = (1-4V2)T            ", "out", [[1,1],[2,-4]]],
 ["Bout  = (1-5V2+4V4)T        ", "out", [[1,1],[2,-5],[4,4]]],
 ["T2    = (1-2V2)T            ", "out", [[1,1],[2,-2]]]
];}

hostrep(nm, H, Av) = {
  print("\n================ HOST ", nm, " ================");
  print("A: ", Av + O(q^9));
  for(i=1,#CL,
    my(cl=CL[i], pl=cl[3], base=if(cl[2]=="in",Sin,Tout));
    my(P0=Peval(pl,0), P1=Peval(pl,1), P2=Peval(pl,2), P3=Peval(pl,3));
    my(xi = if(cl[2]=="in", -P2*G/2, P2*Z2/2));
    my(Bv = Bof(H, mkPhi(base,pl)));
    print("\n--- ", cl[1], "  P(0)=",P0," P(1)=",P1," P(2)=",P2," P(3)=",P3);
    print("    predicted xi_inf = ", xi);
    print("    b_n/a_n at n=10,20,30,",NT-2,": ",
      [polcoeff(Bv,10)/polcoeff(Av,10)*1., polcoeff(Bv,20)/polcoeff(Av,20)*1.,
       polcoeff(Bv,30)/polcoeff(Av,30)*1., polcoeff(Bv,NT-2)/polcoeff(Av,NT-2)*1.]);
    my(Cv = Bv - xi*Av);
    print("    |c_n|^(1/n) for c=B-xi*A, n=40,60,80,",NT-2,": ",
      [abs(polcoeff(Cv,40))^(1./40), abs(polcoeff(Cv,60))^(1./60),
       abs(polcoeff(Cv,80))^(1./80), abs(polcoeff(Cv,NT-2))^(1./(NT-2))]);
    print("    |b_n|^(1/n)               n=40,60,80,",NT-2,": ",
      [abs(polcoeff(Bv,40))^(1./40), abs(polcoeff(Bv,60))^(1./60),
       abs(polcoeff(Bv,80))^(1./80), abs(polcoeff(Bv,NT-2))^(1./(NT-2))]);
    print("    v2(b_n)/n  n=40,60,80,",NT-2,": ",
      [valuation(polcoeff(Bv,40),2)/40., valuation(polcoeff(Bv,60),2)/60.,
       valuation(polcoeff(Bv,80),2)/80., valuation(polcoeff(Bv,NT-2),2)/(NT-2.)]);
  );
};
