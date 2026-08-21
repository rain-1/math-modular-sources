PRECSET=if(type(PS)=="t_INT",PS,120);
read("setup.gp");
u=tC; v=Vd(tC,2); x=tF; G=(1+3*v)/(1-3*u);
\\ reversion: express series in q as series in x
inx(f)={my(r=serreverse(x)); subst(f,q,r)};
