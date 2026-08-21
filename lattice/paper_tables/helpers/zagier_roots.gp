parZ=[[7,2,-8],[9,3,27],[10,3,9],[11,3,-1],[12,4,32],[17,6,72]];
fam=["A","B","C","D","E","F"];
for(i=1,6, my(a=parZ[i][1],c=parZ[i][3]); print(fam[i]," roots(x^2-",a,"x+",c,")=",polroots(x^2-a*x+c)~));
\q
