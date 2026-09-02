default(parisize,"48G");
read("rec_17_15.gp");
print("(17,15): #PV=",#PV,"  degs=",vector(#PV,i,poldegree(PV[i])),"  max digits=",vecmax(vector(#PV,i,vecmax(vector(poldegree(PV[i])+1,j,#Str(abs(polcoeff(PV[i],j-1))))))));
print("  P_0 valuation at n=0: ",valuation(PV[1],n));
read("rec_23_10.gp");
print("(23,10): #PV=",#PV,"  degs=",vector(#PV,i,poldegree(PV[i])),"  max digits=",vecmax(vector(#PV,i,vecmax(vector(poldegree(PV[i])+1,j,#Str(abs(polcoeff(PV[i],j-1))))))));
print("  P_0 valuation at n=0: ",valuation(PV[1],n));
quit;
