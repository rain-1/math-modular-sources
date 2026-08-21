default(realprecision, 60);
G=znstar(12,1);
mf=mfinit([12,3,[G,znconreylog(G,5)]],0);
EB=mfeigenbasis(mf);
print("#EB=",#EB);
for(e=1,#EB, print("f",e," coefs=",mfcoefs(EB[e],10)));
lf=lfunmf(mf,EB[1]);
print("type lf=",type(lf), " #lf=", if(type(lf)=="t_VEC",#lf,1));
if(type(lf)=="t_VEC", for(u=1,#lf, print("u",u," L(1)=",lfun(lf[u],1)," L(2)=",lfun(lf[u],2)," L(3)=",lfun(lf[u],3))),
  print("L(2)=",lfun(lf,2)));
quit
