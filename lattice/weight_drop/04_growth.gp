default(realprecision,40);
NQ = 200;
read("lattice/weight_drop/03_setup.gp"); build();
{for(r=1,#DAT, my(nm=DAT[r][1],N=DAT[r][2],PHI=DAT[r][5]);
  my(g=vector(6,j, my(m=140+10*j); abs(polcoeff(PHI,m))^(1.0/m)));
  print(nm,"  N=",N,"  |c(m)|^(1/m) at m=150..200: ",g);
  print("   1/|q_0| needed < ", exp(2*Pi/sqrt(N)), "   (=e^{2pi/sqrt N}); ratio at m=200: ",
        g[6]*exp(-2*Pi/sqrt(N)));
);}
quit;
