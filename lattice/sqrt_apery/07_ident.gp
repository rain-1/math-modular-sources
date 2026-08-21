default(parisizemax,6000000000);
default(realprecision,300);
xi=0.10018744922933940616775868213306112162877772050081221171621939873497220871205775808509584097889068810390680630219845478301693292894533577046669203809845145295550789880948045187971132622476722346776293367425062982491759259498666437881840066684105157691347522541434986715960370511832050544706080093849818552693848253858164399682340040926845880060698875406784332360763077897549245980741967370291789059664110760732894025170552820558400452719124599668853130221754715116515578458360137153001511885110170388379113841439891342722055287689388410074349497290933896897;
print("algdep deg 1..10:");
for(d=1,10, my(v=algdep(xi,d)); print("  d=",d,"  ",v,"  height=",vecmax(abs(Vec(v)))));
L(D,s)=lfun(lfuncreate(D),s);
s2=sqrt(2); s3=sqrt(3); s6=sqrt(6); l12=log(1+s2);
G=L(-4,2);
b1=[xi,1,zeta(2),zeta(3),Pi^2,Pi^3,G,L(-3,2),L(-8,2),L(-24,2),L(8,2),L(12,2),L(24,2),Pi*log(2),Pi*log(3),Pi*l12];
print("basis1: ",lindep(b1,120)~);
\\ CM periods (Chowla-Selberg) for discriminants -3,-4,-8,-24
CS(D)=my(h,w,P=1); h=qfbclassno(D); w=if(D==-3,6,D==-4,4,2);
  for(j=1,abs(D)-1, if(gcd(j,abs(D))==1, P*=gamma(j/abs(D))^kronecker(D,j)));
  (P)^(w/(4*h))/sqrt(2*Pi*abs(D));
Om24=CS(-24); Om8=CS(-8); Om3=CS(-3); Om4=CS(-4);
print("Omega_-24=",Om24,"  Omega_-8=",Om8);
{for(j=-3,3, print("  xi/(Pi^",j,"*Om24^4) : ",lindep([xi, Pi^j*Om24^4],120)~));}
{for(j=-3,3, print("  xi/(Pi^",j,"*Om8^4) : ",lindep([xi, Pi^j*Om8^4],120)~));}
b3=[xi,Om24^4,Pi*Om24^4,Om24^4/Pi,Om24^2,Pi^2*Om24^2,1,Pi^2];
print("basis3(CM -24): ",lindep(b3,120)~);
b4=[xi,Om8^4,Pi*Om8^4,Om8^4/Pi,Om3^4,Pi*Om3^4,1];
print("basis4(CM -8,-3): ",lindep(b4,120)~);
\\ gamma-quotient style: xi vs Gamma(k/24) products
print("gamma(1/24)-type: ",lindep([xi, gamma(1/24)*gamma(5/24)*gamma(7/24)*gamma(11/24)/(gamma(13/24)*gamma(17/24)*gamma(19/24)*gamma(23/24)),1,Pi],80)~);
\q
