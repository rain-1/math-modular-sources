default(parisize,"3G"); default(realprecision,50);
myav=[[1,9,4],[1,9,8],[1,9,12],[1,9,16],[1,9,20],[1,9,24],[1,9,48],[1,25,4],[1,25,8],[1,25,12],[4,8,12]];
mynm=["E_1","E_2","E_3","E_4","E_5","E_6","E_12","F_1","F_2","F_3","G"];
mk(av)=my(p=prod(i=1,#av,(1-av[i]*x)),a=polcoeff(p,3),b=polcoeff(p,2),c=polcoeff(p,1)); ellinit([0,b,0,a*c,a^2]);
for(i=1,#myav, my(E0=mk(myav[i]),ch,Em,a); a=polcoeff(prod(k=1,#myav[i],(1-myav[i][k]*x)),3); Em=ellminimalmodel(E0,&ch); print(mynm[i]," | N=",ellglobalred(Em)[1]," | j=",Em.j," | ord(x=0 pt)=",ellorder(Em,ellchangepoint([0,a],ch))," | ap=",vector(12,n,ellap(Em,prime(n)))));
print("");
for(k=1,24, my(F,E); F=ellfromeqn(x^2*y + y + x*y^2 + x + k*x*y); E=ellinit(F); if(type(E)=="t_VEC" && #E>0, E=ellminimalmodel(E); print("Boyd k=",k," | N=",ellglobalred(E)[1]," | j=",E.j," | ap=",vector(12,n,ellap(E,prime(n))))));
quit;
