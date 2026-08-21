default(parisizemax, 6000000000);
read("identlib.gp");
default(realprecision, 360);
x = eval(externstr("cat _lim136.txt")[1]);
B = List();
listput(B,[1,"1"]); listput(B,[zeta(2),"zeta(2)"]); listput(B,[zeta(3),"zeta(3)"]);
listput(B,[zeta(4),"zeta(4)"]);
listput(B,[lfun(-4,2),"G"]); listput(B,[lfun(-3,2),"L(2,-3)"]); listput(B,[lfun(-3,3),"L(3,-3)"]);
listput(B,[lfun(-4,3),"L(3,-4)"]); listput(B,[lfun(-8,2),"L(2,-8)"]); listput(B,[lfun(-8,3),"L(3,-8)"]);
listput(B,[lfun(8,2),"L(2,8)"]); listput(B,[lfun(8,3),"L(3,8)"]); listput(B,[lfun(-7,2),"L(2,-7)"]);
listput(B,[lfun(5,2),"L(2,5)"]); listput(B,[Pi^3,"Pi^3"]); listput(B,[Pi^2*log(2),"Pi^2log2"]);
listput(B,[log(2)^3,"log2^3"]); listput(B,[sqrt(2),"sqrt2"]); listput(B,[sqrt(2)*zeta(2),"sqrt2 zeta2"]);
listput(B,[sqrt(2)*zeta(3),"sqrt2 zeta3"]); listput(B,[log(1+sqrt(2)),"log(1+sqrt2)"]);
listput(B,[log(1+sqrt(2))^2,"log(1+sqrt2)^2"]); listput(B,[log(1+sqrt(2))^3,"log(1+sqrt2)^3"]);
listput(B,[Catalan,"Catalan"]);
\\ cusp form L-values weight 3 and 4, levels dividing 64 and 32, s=1,2,3
CB = List();
{
for(li=1,8, my(L=[8,16,32,64,4,12,24,48][li]);
  for(k=3,4,
    my(G=znstar(L,1), seen=vectorsmall(L));
    for(m=1,L,
      if(gcd(m,L)!=1 || seen[m], next);
      my(o=znorder(Mod(m,L)));
      for(j=1,o, if(gcd(j,o)==1, seen[lift(Mod(m,L)^j)]=1));
      if(chareval(G,znconreylog(G,m),-1) != if(k%2==0,0,1/2), next);
      if(eulerphi(o)>2, next);
      iferr(
        my(mf=mfinit([L,k,[G,znconreylog(G,m)]],0));
        if(mfdim(mf)>0,
          my(EB=mfeigenbasis(mf));
          for(e=1,#EB,
            my(lf=lfunmf(mf,EB[e]));
            lf = iferr(lfun(lf,2); [lf], E, lf);
            if(type(lf)!="t_VEC", lf=[lf]);
            for(u=1,#lf, for(s=1,k-1,
              my(v=lfun(lf[u],s));
              listput(CB,[real(v),Str("L(f",L,"w",k,"c",m,"n",e,"e",u,",",s,")")])))))
      ,E,))));
}
B = Vec(B); CB = Vec(CB);
print("#B=",#B," #CB=",#CB);
print("singles(const): ", tryrel(x, B, 350));
print("singles(cusp): ", tryrel(x, concat([[1,"1"]],CB), 350));
print("pairs(const): ", trypair(x, B, 350));
print("pairs(cusp x const): ", trypair2(x, CB, B, 350));
quit
