default(parisize, 2000000000);
gettime();
M = mfinit([28,77/2]);
print("dim M_{77/2}(28) = ", mfdim(M), "  time=", gettime());
M2 = mfinit([28,53/2]);
print("dim M_{53/2}(28) = ", mfdim(M2), "  time=", gettime());
quit;
