default(parisizemax, 6000000000);
default(realprecision, 540);
x = eval(externstr("cat _lim136b.txt")[1]);
best = 0;
{
for(li=1,14, my(L=[3,6,9,12,18,24,36,72,144,27,54,108,5,10][li]);
  for(k=3,5,
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
              my(v=real(lfun(lf[u],s)));
              if(v==0, next);
              my(vv=lindep([1,x,v]));
              if(#vv==3 && vv[2]!=0 && vecmax(abs(vv))<100000 && abs(vv[1]+vv[2]*x+vv[3]*v)<10.0^(-500),
                 print("MATCH L(f_",L,"w",k,"c",m,"n",e,"e",u,",",s,") : ",vv); best=1)))))
      ,E,))));
}
print("done, found=",best);
quit
