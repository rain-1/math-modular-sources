\\ C1-C5, C7, C12 : integral monodromy checks
mz0 = [-1,2,-1,-1; 1,-4,2,2; 1,8,-4,-3; -6,-19,9,7];
mzi = [-8,-5,2,1; 14,9,-4,-2; -22,-15,5,3; 45,32,-12,-7];
mz1 = [2,1,0,0; -3,-2,0,0; 0,0,1,0; 0,0,-3,-1];
jz  = [0,0,0,1; 0,0,1,1; 0,-1,0,0; -1,-1,0,0];
id4 = matid(4);

print("=== C1 symplectic ===");
print("det jz = ", matdet(jz));
print("jz alternating? jz+jz~ = ", jz+mattranspose(jz));
{
foreach([["mz0",mz0],["mzi",mzi],["mz1",mz1]], p,
  my(nm=p[1], mm=p[2]);
  print(nm, ": det = ", matdet(mm), "   mm^T*jz*mm - jz = ", mattranspose(mm)*jz*mm - jz);
);
}

print();
print("=== C2 char polys / Jordan ===");
{
foreach([["mz0",mz0],["mzi",mzi],["mz1",mz1]], p,
  my(nm=p[1], mm=p[2]);
  print(nm, " charpoly = ", factor(charpoly(mm,'x)));
  print(nm, " minpoly  = ", factor(minpoly(mm,'x)));
);
}
print("rank(mz0^3-I) = ", matrank(mz0^3-id4), "   mz0^3-I = ", mz0^3-id4);
print("mz0 order? mz0^3 = ", mz0^3);
print("mz1^2 - I = ", mz1^2-id4);
print("rank(mzi^6-I) = ", matrank(mzi^6-id4), "   mzi^6-I = ", mzi^6-id4);
print("rank(mzi^2-I) = ", matrank(mzi^2-id4));
print("(mzi+I) rank = ", matrank(mzi+id4), " ; (mzi+I)^2 rank = ", matrank((mzi+id4)^2));
print("(mz0^3-I)^2 = ", (mz0^3-id4)^2);

print();
print("=== C3 product relations ===");
{
my(L=[["mz0",mz0],["mz1",mz1],["mzi",mzi]]);
for(a=1,3, for(b=1,3, for(c=1,3,
  if(a!=b && b!=c && a!=c,
    for(e1=0,1, for(e2=0,1, for(e3=0,1,
      my(A=if(e1,L[a][2]^-1,L[a][2]), B=if(e2,L[b][2]^-1,L[b][2]), C=if(e3,L[c][2]^-1,L[c][2]));
      my(P=A*B*C);
      if(P==id4,
        print("IDENTITY: ",L[a][1],if(e1,"^-1",""),"*",L[b][1],if(e2,"^-1",""),"*",L[c][1],if(e3,"^-1",""))
      );
      if(P==-id4,
        print("MINUS-IDENTITY: ",L[a][1],if(e1,"^-1",""),"*",L[b][1],if(e2,"^-1",""),"*",L[c][1],if(e3,"^-1",""))
      );
    )));
  );
)));
}
print("-- traces of all 6 unsigned products (no inverses):");
{
my(L=[["mz0",mz0],["mz1",mz1],["mzi",mzi]]);
for(a=1,3, for(b=1,3, for(c=1,3,
  if(a!=b && b!=c && a!=c,
    my(P=L[a][2]*L[b][2]*L[c][2]);
    print(L[a][1],"*",L[b][1],"*",L[c][1]," = ",P, "  charpoly=",charpoly(P,'x));
  );
)));
}
print("-- is any product conjugate-to-identity? check orders:");
{
my(L=[["mz0",mz0],["mz1",mz1],["mzi",mzi]]);
for(a=1,3, for(b=1,3, for(c=1,3,
  if(a!=b && b!=c && a!=c,
    my(P=L[a][2]*L[b][2]*L[c][2], o=0, Q=P);
    for(k=1,60, if(Q==id4, o=k; break); Q=Q*P);
    print(L[a][1],"*",L[b][1],"*",L[c][1]," order=",if(o,o,"infinite/>60"));
  );
)));
}

print();
print("=== C4 ===");
n0 = mz0^3 - id4;
print("N0 = ", n0);
print("N0^2 = ", n0^2, "  rank N0 = ", matrank(n0));
s0a = jz*n0; s0b = mattranspose(n0)*jz;
print("S0a = jz*N0 = ", s0a, "  symmetric? ", s0a==mattranspose(s0a));
print("S0b = N0^T*jz = ", s0b, "  symmetric? ", s0b==mattranspose(s0b));
print("S0a eigen (charpoly) = ", factor(charpoly(s0a,'x)));
print("SNF(S0a) = ", matsnf(s0a), "  SNF(S0b) = ", matsnf(s0b));
print("rank S0a = ", matrank(s0a));
print("S0a + S0b = ", s0a+s0b);
quit();
