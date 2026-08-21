\\ Period-annihilation (Dirichlet-polynomial) facts from the book. Run: gp -q verify/period_annihilation.gp
P(c,d,s)=sum(i=1,#c,c[i]*d[i]^(-s));
chk(name,cond)=print(if(cond,"PASS ","FAIL "),name);
\\ Beukers zeta(3), level 6
c=[1,-28,63,-36];d=[1,2,3,6]; chk("Beukers P(0)=P(2)=P(4)=0",vector(3,j,P(c,d,2*j-2))==[0,0,0]);
\\ level-8 zeta(3) full projector = (1-V2)(1-4V2)(1-16V2)
chk("level-8 z3 = prod(1-4^a X)",Pol([-64,84,-21,1])==prod(a=0,2,(1-4^a*x)));
\\ zeta(5) level-16
c=[1,-85,1428,-5440,4096];d=[1,2,4,8,16];
chk("z5 level16 kills s=0,2,4,6",vector(4,j,P(c,d,2*j-2))==[0,0,0,0]);
chk("z5 Gaussian binomial form",Polrev(c)==prod(a=0,3,(1-4^a*x)));
\\ zeta(7) level-12 anti-Fricke parent, target 209/1728
c=[1,-572,11583,-36608,46332,-20736];d=[1,2,3,4,6,12];
chk("z7 level12 kills s=0,2,4,6,8",vector(5,j,P(c,d,2*j-2))==[0,0,0,0,0]);
chk("z7 target -P(7)/2 = 209/1728",-P(c,d,7)/2==209/1728);
\\ beta(4) level 24: inner = D C_B D^-1
c=[1,-56,189,-216];d=[1,2,3,6];
chk("beta4 inner kills s=1,3,5",vector(3,j,P(c,d,2*j-1))==[0,0,0]);
chk("beta4 inner P(4)=-1/3",P(c,d,4)==-1/3);
chk("beta4 inner = Mellin shift of Beukers",c==vector(4,i,[1,-28,63,-36][i]*d[i]));
chk("beta4 level-16 Vandermonde det = -135/8192",matdet([1,1/2,1/4;1,1/8,1/64;1,1/32,1/1024])==-135/8192);
\\ L(2,chi_-3) benchmark: (1-2V2)(1-8V2) = (1,-10,16)
chk("chi-3 purified (1,-10,16)",Vec((1-2*x)*(1-8*x))==[16,-10,1]);
