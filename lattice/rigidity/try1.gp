PRECSET=70;
read("setup.gp"); read("ratfit.gp");
u=tC; v=Vd(tC,2);
scan(u,tF,"tC(tau) in tF",4,55);
scan(v,tF,"tC(2tau) in tF",4,55);
scan(FF/FC,tF,"FF/FC in tF",4,55);
