read("common.gp");
gettime();
BB=build(3,7,3,81,700);
print("delta build NT=700: ",gettime()," ms");
print("t coeffs at 100,300,500,700: ", [polcoeff(truncate(BB[1]),100),polcoeff(truncate(BB[1]),300),polcoeff(truncate(BB[1]),500),polcoeff(truncate(BB[1]),700)]);
print("F coeffs at 100,700: ", [polcoeff(truncate(BB[2]),100),polcoeff(truncate(BB[2]),700)]);
BB=build(3,11,5,125,700); print("eta build: ",gettime()," ms  t[700]=",polcoeff(truncate(BB[1]),700));
BB=build(2,9,3,27,700); print("B build: ",gettime()," ms  t[700]=",polcoeff(truncate(BB[1]),700));
quit;
