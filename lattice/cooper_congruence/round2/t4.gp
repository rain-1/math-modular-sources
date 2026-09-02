read("wt2b.gp");
default(realprecision,50);
t = 0.11+0.29*I;
print("fhat via q-series : ", fhatq(t,7,FS7,300));
print("(compare with the E2*/E4 route: 2.5439318448907179159 + 10.027618844160864517*I)");
print("at the elliptic point (-5+I*sqrt(3))/14 : ", fhatq((-5+I*sqrt(3))/14,7,FS7,300));
quit;
