\\ 76_Gmod.gp -- mod-p structure of the extended kernel G(d) = c(d)/d (row s7).
L = read("/home/ubuntu/code/math-modular-sources/lattice/cooper_congruence/round2/73_cd_s7_2500.txt");
G = List();
{ for(i=1,#L, if(type(L[i][2])=="t_INT", listput(G,[L[i][1], L[i][2]/L[i][1]]))); }
V = Vec(G);
print("entries: ", #V);
print("G(d) mod 2 for the first 30 d: ");
print(vector(30,i,[V[i][1], V[i][2]%2]));
\\ periodicity of G mod p in d
{ forprime(p=2,7,
    my(found=0);
    for(T=1,600,
      my(ok=1);
      for(i=1,#V, my(j=i+1);
        while(j<=#V && (V[j][1]-V[i][1])%T!=0, j++);
        if(j<=#V, if((V[i][2]-V[j][2])%p!=0, ok=0; break)));
      if(ok, print("  G mod ",p," is periodic in d with period ",T); found=1; break));
    if(!found, print("  G mod ",p,": no period <= 600")));
}
\\ zero set and growth
{ my(z=List()); for(i=1,#V, if(V[i][2]==0, listput(z,V[i][1])));
  print("");
  print("G(d) = 0 exactly at d = ", Vec(z));
  print("  all divisible by 7? ", if(#select(x->x%7!=0, Vec(z))==0, "YES", "no"));
  print("  every admissible d divisible by 7 is in the zero set? ",
        if(#select(x->x[1]%7==0 && x[2]!=0, V)==0, "YES", "no")); }
quit;
