default(realprecision, 200);
read("theta_gp.txt");
xi1 = 0.06786721954287826166646998964563721793157057974785253811627679599873293240387754600472436904131064158841277383301338254691816105592800291978569842592232280720350666634075477524078495373724845256233523035364407420441444991361676976050336957323807362044 - 0.06907819963322225514670799790743491919071785808412459538740013509317488539578332659339145894820192577647267658296914645712671900760850497706419803607235071864675633771938344206632717540610599368109439301877396417892631629513247150535996752102799884096*I;
xi2 = 0.1642087341590574584714225970206947065830122999536287405694316515422586105676807745518762686902555647720236146720735800852953327485050302342736826060045161009106544716484284489377941968395761197077123295366087484842016287741923684399374479825889396610 - 0.2369659751659112010711143299806665676669665318863734518080024379160442279419310559148346593990524891514289658806059918517094079894759328173258155337192594168909099387369602162580142247035701450032542277433112273060653991607390105577213119189860336177*I;
g3 = gamma(1/3); g4 = gamma(1/4);
A7 = gamma(1/7)*gamma(2/7)*gamma(4/7)/(gamma(3/7)*gamma(5/7)*gamma(6/7));
alghit(zz, dmax, ht) = { for(dd=1,dmax, my(pp=algdep(zz,dd)); if(poldegree(pp)>0 && vecmax(abs(Vec(pp))) < ht, my(rt=polroots(pp)); my(mm=vecmin(vector(#rt,ii,abs(rt[ii]-zz)))); if(mm < 1e-110, return([dd,pp,mm])))); 0; }
rat(tag, zz) = { my(h=alghit(zz, 26, 10^8)); if(h!=0, print("  RATIO-ALG ", tag, " deg ",h[1]," ",h[2]," resid ",h[3]), print("  ratio ", tag, ": no algebraic relation deg<=26 ht<1e8")); }
rat("Re(xi1)/Im(xi1)", real(xi1)/imag(xi1));
rat("Re(xi2)/Im(xi2)", real(xi2)/imag(xi2));
rat("Re(th1)/Im(th1)", real(th1)/imag(th1));
rat("Re(th2)/Im(th2)", real(th2)/imag(th2));
rat("Re(xi1)/Re(th1)", real(xi1)/real(th1));
rat("Im(xi1)/Im(th1)", imag(xi1)/imag(th1));
rat("Re(xi2)/Re(th2)", real(xi2)/real(th2));
rat("Im(xi2)/Im(th2)", imag(xi2)/imag(th2));
sweep(tag, val, gg, gn, den, amin, amax, dmax) = { my(cnt=0);
  for(aa=amin,amax, for(bb=-8,8,
    my(mo = gg^(aa/den) * Pi^bb); my(zz = val/mo); my(h = alghit(zz, dmax, 10^7));
    if(h != 0, print("  ALG-HIT ", tag, " / (",gn,"^",aa,"/",den," Pi^",bb,") : deg ",h[1],"  minpoly ",h[2],"  resid ",h[3]); cnt++);
  ));
  print(tag, " [",gn,"]: ", cnt, " hits"); }
sweep("Re(xi1)", real(xi1), g3, "G(1/3)", 1, -18, 18, 18);
sweep("Im(xi1)", imag(xi1), g3, "G(1/3)", 1, -18, 18, 18);
sweep("Re(xi1)", real(xi1), A7, "A7", 2, -8, 8, 18);
sweep("Im(xi1)", imag(xi1), A7, "A7", 2, -8, 8, 18);
sweep("Re(th1)", real(th1), g3, "G(1/3)", 1, -18, 18, 18);
sweep("Im(th1)", imag(th1), g3, "G(1/3)", 1, -18, 18, 18);
sweep("Re(xi2)", real(xi2), g4, "G(1/4)", 1, -20, 20, 18);
sweep("Re(th2)", real(th2), g4, "G(1/4)", 1, -20, 20, 18);
sweep("Im(th2)", imag(th2), g4, "G(1/4)", 1, -20, 20, 18);
