default(parisizemax, 6000000000);
NQ = 34;
NA = 30;
useries(dv, r, nq) = q*prod(t=1, #dv, eta(q^dv[t] + O(q^nq))^r[t]);
Fseries(dv, r, nq) = 1 - sum(t=1, #dv, r[t]*dv[t]*sum(n=1, (nq-1)\dv[t], sigma(n)*q^(dv[t]*n))) + O(q^nq);
peel2(Fs, xs, na, nq) = my(a=vector(na+1), G=Fs, xp=1+O(q^nq)); for(n=0, na, my(c=polcoeff(G,n)); a[n+1]=c; if(c!=0, G=G-c*xp); xp=xp*xs); a;
fitrecR(av, ord, dg) = my(nv=(ord+1)*(dg+1), rows=List()); for(n=0, #av-ord-2, my(row=vector(nv)); for(j=0,ord, for(e=0,dg, row[j*(dg+1)+e+1] = n^e*av[n+j+1])); listput(rows,row)); matker(matconcat(Vec(rows)~));

dv=[1, 2]; r=[-24, 24]; C=4096;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(2,"|",4096,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3]; r=[-12, 12]; C=729;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(3,"|",729,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4]; r=[-8, 0, 8]; C=256;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(4,"|",256,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 5]; r=[-6, 6]; C=125;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(5,"|",125,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6]; r=[-4, -4, 4, 4]; C=81;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(6,"|",81,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6]; r=[-5, 1, -1, 5]; C=72;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(6,"|",72,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6]; r=[-6, 6, -6, 6]; C=64;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(6,"|",64,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 7]; r=[-4, 4]; C=49;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(7,"|",49,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[8, -40, 40, -8]; C=256;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",256,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[4, -26, 26, -4]; C=128;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",128,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[0, -12, 12, 0]; C=64;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",64,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[-4, 2, -2, 4]; C=32;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",32,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[-8, 16, -16, 8]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",16,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[-12, 30, -30, 12]; C=8;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",8,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8]; r=[-16, 44, -44, 16]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(8,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 9]; r=[-3, 0, 3]; C=27;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(9,"|",27,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 5, 10]; r=[-2, -2, 2, 2]; C=25;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(10,"|",25,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 5, 10]; r=[-3, 1, -1, 3]; C=20;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(10,"|",20,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 5, 10]; r=[-4, 4, -4, 4]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(10,"|",16,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[3, -12, -9, 9, 12, -3]; C=64;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",64,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[6, -20, -10, 10, 20, -6]; C=144;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",144,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[2, -10, -6, 6, 10, -2]; C=48;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",48,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[-2, 0, -2, 2, 0, 2]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",16,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[9, -28, -11, 11, 28, -9]; C=324;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",324,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[5, -18, -7, 7, 18, -5]; C=108;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",108,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[1, -8, -3, 3, 8, -1]; C=36;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",36,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[-3, 2, 1, -1, -2, 3]; C=12;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",12,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[-7, 12, 5, -5, -12, 7]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",4,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[4, -16, -4, 4, 16, -4]; C=81;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",81,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[0, -6, 0, 0, 6, 0]; C=27;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",27,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[-4, 4, 4, -4, -4, 4]; C=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",9,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[-8, 14, 8, -8, -14, 8]; C=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",3,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 12]; r=[-12, 24, 12, -12, -24, 12]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(12,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 13]; r=[-2, 2]; C=13;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(13,"|",13,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 7, 14]; r=[-3, 3, -3, 3]; C=8;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(14,"|",8,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 5, 15]; r=[-2, 2, -2, 2]; C=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(15,"|",9,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[4, -14, 0, 14, -4]; C=64;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",64,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[2, -9, 0, 9, -2]; C=32;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",32,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[0, -4, 0, 4, 0]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",16,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[-2, 1, 0, -1, 2]; C=8;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",8,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[-4, 6, 0, -6, 4]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",4,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[-6, 11, 0, -11, 6]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",2,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16]; r=[-8, 16, 0, -16, 8]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(16,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[2, -4, -10, 10, 4, -2]; C=36;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",36,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[1, -2, -9, 9, 2, -1]; C=24;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",24,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[0, 0, -8, 8, 0, 0]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",16,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[-1, -1, 0, 0, 1, 1]; C=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",9,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[-2, 1, 1, -1, -1, 2]; C=6;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",6,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[-3, 3, 2, -2, -3, 3]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",4,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 9, 18]; r=[-6, 6, 12, -12, -6, 6]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(18,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 10, 20]; r=[2, -8, 2, -2, 8, -2]; C=25;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(20,"|",25,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 10, 20]; r=[-2, 2, -2, 2, -2, 2]; C=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(20,"|",5,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 10, 20]; r=[-6, 12, -6, 6, -12, 6]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(20,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 10, 20]; r=[1, -6, 5, -5, 6, -1]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(20,"|",16,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 7, 21]; r=[-1, -1, 1, 1]; C=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(21,"|",7,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 11, 22]; r=[-2, 2, -2, 2]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(22,"|",4,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[0, 1, -4, -7, 7, 4, -1, 0]; C=12;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",12,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[1, -2, -3, -6, 6, 3, 2, -1]; C=18;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",18,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-2, 4, -2, -4, 4, 2, -4, 2]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-1, 1, -1, -3, 3, 1, -1, 1]; C=6;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",6,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[0, -2, 0, -2, 2, 0, 2, 0]; C=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",9,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[1, -3, -3, -1, 1, 3, 3, -1]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",16,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[2, -6, -2, 0, 0, 2, 6, -2]; C=24;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",24,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[3, -9, -1, 1, -1, 1, 9, -3]; C=36;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",36,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-3, 4, 1, 0, 0, -1, -4, 3]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",2,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-2, 1, 2, 1, -1, -2, -1, 2]; C=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",3,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[0, -3, 0, 3, -3, 0, 3, 0]; C=8;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",8,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[1, -6, 1, 4, -4, -1, 6, -1]; C=12;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",12,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-4, 4, 4, 4, -4, -4, -4, 4]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 24]; r=[-1, -3, 3, 7, -7, -3, 3, 1]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(24,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 5, 25]; r=[-1, 0, 1]; C=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(25,"|",5,"|",B,"|",1,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 9, 27]; r=[0, -4, 4, 0]; C=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(27,"|",9,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 7, 14, 28]; r=[0, -2, 0, 0, 2, 0]; C=7;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(28,"|",7,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 7, 14, 28]; r=[-4, 8, -4, 4, -8, 4]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(28,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 7, 14, 28]; r=[-1, 0, 1, -1, 0, 1]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(28,"|",4,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 5, 6, 10, 15, 30]; r=[-1, 1, -1, -1, 1, 1, -1, 1]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(30,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 5, 6, 10, 15, 30]; r=[-2, 2, 1, 1, -1, -1, -2, 2]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(30,"|",2,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 5, 6, 10, 15, 30]; r=[-3, 3, 3, 3, -3, -3, -3, 3]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(30,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16, 32]; r=[0, 0, -6, 6, 0, 0]; C=8;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(32,"|",8,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16, 32]; r=[2, -7, 3, -3, 7, -2]; C=16;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(32,"|",16,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16, 32]; r=[-2, 3, -1, 1, -3, 2]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(32,"|",2,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16, 32]; r=[0, -4, 8, -8, 4, 0]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(32,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 11, 33]; r=[-1, 1, -1, 1]; C=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(33,"|",3,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[1, -4, 0, 1, 0, -1, 0, 4, -1]; C=9;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",9,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[-1, 1, 0, -1, 0, 1, 0, -1, 1]; C=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",3,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[-3, 6, 0, -3, 0, 3, 0, -6, 3]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[1, -2, -3, 0, 0, 0, 3, 2, -1]; C=12;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",12,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[-1, 3, -3, -2, 0, 2, 3, -3, 1]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[-2, 0, 4, 2, 0, -2, -4, 0, 2]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 9, 12, 18, 36]; r=[0, -3, 1, 3, 0, -3, -1, 3, 0]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(36,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 8, 10, 20, 40]; r=[0, -2, 2, 0, 0, -2, 2, 0]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(40,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 8, 10, 20, 40]; r=[1, -3, -1, -1, 1, 1, 3, -1]; C=10;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(40,"|",10,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 8, 10, 20, 40]; r=[-1, 0, 2, 1, -1, -2, 0, 1]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(40,"|",2,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 8, 10, 20, 40]; r=[0, -1, -1, 0, 0, 1, 1, 0]; C=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(40,"|",5,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 8, 10, 20, 40]; r=[-2, 2, 2, 2, -2, -2, -2, 2]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(40,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 7, 14, 21, 42]; r=[-2, 2, 2, -2, 2, -2, -2, 2]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(42,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 5, 9, 15, 45]; r=[0, -2, 0, 0, 2, 0]; C=5;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(45,"|",5,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 23, 46]; r=[-1, 1, -1, 1]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(46,"|",2,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]; r=[0, -1, 0, 0, -1, 1, 0, 0, 1, 0]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(48,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]; r=[-1, 1, 1, -1, -2, 2, 1, -1, -1, 1]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(48,"|",2,"|",B,"|",2,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]; r=[-2, 3, 2, -2, -3, 3, 2, -2, -3, 2]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(48,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]; r=[1, -4, -1, 3, 3, -3, -3, 1, 4, -1]; C=6;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(48,"|",6,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 12, 16, 24, 48]; r=[0, -2, 0, 2, 2, -2, -2, 0, 2, 0]; C=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(48,"|",3,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 13, 26, 52]; r=[-2, 4, -2, 2, -4, 2]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(52,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 7, 8, 14, 28, 56]; r=[-1, 2, -2, -1, 1, 2, -2, 1]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(56,"|",2,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 5, 6, 10, 12, 15, 20, 30, 60]; r=[-1, 0, 1, 1, 1, 0, 0, -1, -1, -1, 0, 1]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(60,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 7, 9, 21, 63]; r=[-1, 2, 1, -1, -2, 1]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(63,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16, 32, 64]; r=[0, 0, -2, 0, 2, 0, 0]; C=4;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(64,"|",4,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 16, 32, 64]; r=[0, -2, 3, 0, -3, 2, 0]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(64,"|",2,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 5, 7, 10, 14, 35, 70]; r=[-1, 1, 1, 1, -1, -1, -1, 1]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(70,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 9, 12, 18, 24, 36, 72]; r=[-1, 1, 0, 1, 0, -1, 1, 0, -1, 0, -1, 1]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(72,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 4, 6, 8, 9, 12, 18, 24, 36, 72]; r=[-1, 2, 2, -2, -6, 1, -1, 6, 2, -2, -2, 1]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(72,"|",2,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 3, 6, 13, 26, 39, 78]; r=[-1, 1, 1, -1, 1, -1, -1, 1]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(78,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 3, 9, 27, 81]; r=[0, -1, 0, 1, 0]; C=3;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(81,"|",3,"|",B,"|",3,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 8, 11, 22, 44, 88]; r=[0, -1, 1, 0, 0, -1, 1, 0]; C=2;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(88,"|",2,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

dv=[1, 2, 4, 5, 10, 20, 25, 50, 100]; r=[-1, 2, -1, 0, 0, 0, 1, -2, 1]; C=1;
us = useries(dv,r,NQ); Fs = Fseries(dv,r,NQ);
for(B=1, 90, my(xs = us/(1+B*us+C*us^2)); my(a = peel2(Fs, xs, NA, NQ)); my(ker=fitrecR(a,2,3)); if(matsize(ker)[2]==1, my(P=vector(3,j,sum(e=0,3,ker[(j-1)*4+e+1,1]*x^e))); my(cp=sum(j=1,3,polcoeff(P[j],3)*y^(j-1))); my(lead=polcoeff(cp,2)); if(lead!=0, my(cc=polcoeff(cp,0)/lead, ss=-polcoeff(cp,1)/lead); print(100,"|",1,"|",B,"|",4,"|",ss,"|",cc,"|",matsize(ker)[2]))));

quit;
