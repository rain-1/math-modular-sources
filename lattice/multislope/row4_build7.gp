default(parisizemax, 6000000000);
default(realprecision,120);
NORD = 700;

fk(k,NN) = { my(r=1+O(q^(NN+1))); for(n=1, NN\k+1, r*=(1-q^(k*n))); r }

rq = q*(fk(2,NORD)*fk(12,NORD)/(fk(4,NORD)*fk(6,NORD)))^6;
zq = rq/(1+rq);
Qz = serreverse(zq);
print("serreverse done");

cd = [1,-588,11583,-27456,-138996,585728,-762048,331776];
dl = [1,2,3,4,6,8,12,24];
gvec = vector(NORD, m, sum(i=1,#dl, if(m%dl[i]==0, cd[i]*sigma(m\dl[i],7), 0)));
PSIq = sum(m=1,NORD, (gvec[m]/m^7)*q^m) + O(q^(NORD+1));

Ms = NORD\2+3;
apery = vector(Ms+1, mp1, my(m=mp1-1); sum(k=0,m, binomial(m,k)^2*binomial(m+k,k)^2));
Aser_s = sum(m=0,Ms, apery[m+1]*'s^m) + O('s^(Ms+1));
sqrtpart = sqrt(1 - 34*'s + 's^2 + O('s^(Ms+1)));
Acubed_s = Aser_s^3 * sqrtpart;
print("A(s) done");

zser = q + O(q^(NORD+1));
sser = (zser/(1-zser))^2;
Aser = subst(Acubed_s, 's, sser);
PSIofQz = subst(PSIq, q, Qz);
Bser = Aser * PSIofQz;
print("composition done");

NN = NORD-3;
An = vector(NN+1, i, polcoeff(Aser, i-1));
Bn = vector(NN+1, i, polcoeff(Bser, i-1));

print("NN = ", NN);
print("A_0..A_12 = ", vector(13,i,An[i]));
print("B_1..B_8  = ", vector(8,i,Bn[i+1]));
{my(ok=1); for(i=1,NN+1, if(type(An[i])!="t_INT", ok=0; print("A_",i-1," NOT integer")));
 print("all A_n integral? ", ok);}
{my(ok=1, dn=1, kk=0); for(n=1,NN, dn=lcm(dn,n); my(v=dn^7*Bn[n+1]); if(denominator(v)!=1, ok=0; if(kk<5,print("FAIL n=",n); kk++)));
 print("d_n^7 B_n integral for n<=",NN,"? ", ok);}
{my(ok=1, dn=1); for(n=1,NN, dn=lcm(dn,n); my(v=dn^6*Bn[n+1]); if(denominator(v)!=1, ok=0; break));
 print("d_n^6 B_n integral (i.e. k=7 not sharp)? ", ok);}
print("B_n/A_n at n=",NN," = ", 1.0*Bn[NN+1]/An[NN+1]);
print("target = ", 1.0*(1463/13824)*zeta(7));
print("diff = ", 1.0*Bn[NN+1]/An[NN+1] - (1463/13824)*zeta(7));

write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_An7.txt", An);
write("/home/ubuntu/code/math-modular-sources/lattice/multislope/row4_Bn7.txt", Bn);
print("written");
quit
