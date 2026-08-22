default(parisizemax, 4000000000);
sqfp(dd) = my(f=factor(dd), r=1); for(i=1,#f~, if(f[i,2]%2==1, r=r*f[i,1])); if(dd<0 && r>0, r=-r); r
nm = ["1 Zagier B","2 Zagier D","3 Zagier E","4 Herf#4","5 Herf#6","6 sqrt(AZ(7,3,81))","7 Beukers/sqrtApery","8 Apery zeta3","9 AZ eta","10 AZ delta","11 AZ(9,3,-27)","13 AESZ184"];
cps = [x^2-9*x+27, x^2-11*x-1, x^2-12*x+32, x^2-117*x+3969, x^2-72*x+1728, x^2-56*x+1296, x^2-136*x+16, x^2-34*x+1, x^2-22*x+125, x^2-14*x+81, x^2-18*x-27, x^2-88*x+2000];
for(i=1,#cps, my(dd=poldisc(cps[i])); print(nm[i], "  disc=", dd, " = ", if(dd!=0,factor(dd),0), "  root field Q(sqrt(", sqfp(dd), "))  places: ", if(dd>0, if(issquare(dd),"splits over Q (2 real, K=Q)","2 real"), "1 complex")));
print("");
print("Beukers quartic vs Q(zeta_20)^+ : ", polsubcyclo(20,4));
print("  nfisisom ? ", nfisisom(x^4-496*x^3+1006*x^2-496*x+1, polsubcyclo(20,4)));
