default(parisizemax, 6000000000);
read("sc_rows.gp");
{
gettime();
A5 = genseq(R5cf, 4, [1], 400);
print("R5 A N=400: ", gettime(), " ms; digits(A_400)=", #digits(A5[401]));
X51 = compan(R5cf, 4, 1, 400);
print("R5 X1 N=400: ", gettime(), " ms; den digits=", #digits(denominator(X51[401])));
print("check B: ", vector(4,i,X51[i+1]));
A4 = genseq(R4cf, 6, [1], 400);
print("R4 A N=400: ", gettime(), " ms");
X41 = compan(R4cf, 6, 1, 400);
print("R4 X1 N=400: ", gettime(), " ms");
print("R4 X1 first 8: ", vector(8,i,X41[i]));
}
