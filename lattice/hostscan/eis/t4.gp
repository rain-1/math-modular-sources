plist(f) = my(G,r); G=znstar(f,1); r=List(); for(n=1, f, if(gcd(n,f)==1, if(type(znconreyconductor(G,n))=="t_INT", listput(r,[n, zncharisodd(G,n), zncharorder(G,n)])))); Vec(r)
print("f=1 ",plist(1));
print("f=3 ",plist(3));
print("f=4 ",plist(4));
print("f=5 ",plist(5));
print("f=6 ",plist(6));
print("f=8 ",plist(8));
print("f=12 ",plist(12));
G3=znstar(3,1); G12=znstar(12,1);
print("induce ",zncharinduce(G3,2,12));
print("mul ",zncharmul(G12, zncharinduce(G3,2,12), 1));
print("isodd ", zncharisodd(G12, zncharinduce(G3,2,12)));
print("dim7 ", mfdim([7,4,0],3));
print("dim12 ", mfdim([12,4,0],3));
print("dim60 ", mfdim([60,3,0],3));
quit;
