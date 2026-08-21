default(parisizemax,6000000000); default(realprecision,260);
xi=0.10018744922933940616775868213306112162877772050081221171621939873497220871205775808509584097889068810390680630219845478301693292894533577046669203809845145295550789880948045187971132622476722346776293367425062982491759259498666437881840066684105157691347522541434986715960370511832050544706080093849818552693848253858164399682340040926845880060698875406784332360763077897549245980741967370291789059664110760732894025170552820558400452719124599668853130221754715116515578458360137153001511885110170388379113841439891342722055287689388410074349497;
{CS(D)=my(h=qfbclassno(D), w=if(D==-3,6,if(D==-4,4,2)), P=1.0);
 for(j=1,abs(D)-1, if(gcd(j,abs(D))==1, P=P*gamma(j/abs(D))^kronecker(D,j)));
 P^(w/(4*h))/sqrt(2*Pi*abs(D));}
L(D,s)=lfun(lfuncreate(D),s);
S=List(); nm=List();
ad(x,s)=listput(S,x);listput(nm,s);
ad(1,"1"); ad(zeta(3),"zeta3"); ad(zeta(2),"zeta2"); ad(Pi^3,"Pi^3"); ad(Pi^2,"Pi^2");
ad(L(-4,2),"G"); ad(L(-3,2),"L(-3,2)"); ad(L(-8,2),"L(-8,2)"); ad(L(-24,2),"L(-24,2)");
ad(L(8,2),"L(8,2)"); ad(L(12,2),"L(12,2)"); ad(L(24,2),"L(24,2)");
ad(L(-3,3),"L(-3,3)"); ad(L(-4,3),"L(-4,3)"); ad(L(-8,3),"L(-8,3)"); ad(L(-24,3),"L(-24,3)");
ad(Pi*log(2),"Pi log2"); ad(Pi*log(3),"Pi log3"); ad(Pi*log(1+sqrt(2)),"Pi log(1+r2)");
ad(Pi^2*log(2),"Pi^2 log2"); ad(Pi^2*log(1+sqrt(2)),"Pi^2 log(1+r2)");
ad(log(1+sqrt(2))^3,"log(1+r2)^3"); ad(log(2)^3,"log2^3");
{for(k=1,4, my(D=[-3,-4,-8,-24][k], O=CS(D));
  for(j=-2,3, ad(Pi^j*O^4, concat(["Pi^",Str(j)," Om(",Str(D),")^4"]));
              ad(Pi^j*O^2, concat(["Pi^",Str(j)," Om(",Str(D),")^2"]));
              ad(Pi^j*O^6, concat(["Pi^",Str(j)," Om(",Str(D),")^6"]))));}
print("basis size ",#S);
hits=0;
{for(i=1,#S, my(v=lindep([xi,S[i]],200));
  if(v!=0 && vecmax(abs(v))<10^14, hits++; print("  2-term HIT: xi vs ",nm[i],"  ",v~)));}
print("2-term hits: ",hits);
\\ 3-term with 1
hits=0;
{for(i=1,#S, my(v=lindep([xi,S[i],1],160));
  if(v!=0 && vecmax(abs(v))<10^9 && v[1]!=0, hits++; print("  3-term HIT: xi, ",nm[i],", 1 : ",v~)));}
print("3-term hits: ",hits);
\\ sqrt multiples of xi
{for(k=1,6, my(r=[sqrt(2),sqrt(3),sqrt(6),1/sqrt(2),1/sqrt(3),1/sqrt(6)][k]);
  for(i=1,#S, my(v=lindep([xi*r,S[i]],200));
    if(v!=0 && vecmax(abs(v))<10^14, print("  scaled HIT: xi*",r," vs ",nm[i],"  ",v~))));}
print("done");
\q
