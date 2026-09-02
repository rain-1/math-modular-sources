\\ 56_cd_s7chk.gp -- SANITY: reproduce the parallel session's c(d) for s7 at non-square d.
read("56_cdlib.gp");
default(realprecision, 60);
initfser(1,2200);
TGT = [[8,16],[29,145],[37,148],[53,795],[65,650],[109,9919],[113,9266],[1,1],[4,-8],[9,9],[16,32],[25,-125],[21,0],[57,798],[85,2975],[93,2325]];
{
for(j=1,#TGT,
  my(d=TGT[j][1], want=TGT[j][2], bt=chosenbeta(1,d), R, v);
  if(bt<0, print("  d=",d,"  NOT ADMISSIBLE"); next);
  R = Traw(1,d,bt,40);
  if(type(R)=="t_STR", print("  d=",d,"  NOREP"); next);
  v = I*sqrt(3)*R[1];
  print("  d=",d,"  beta=",bt,"  cls=",R[2],"/",R[3],"   c(d)=",v,
        "   want ",want,"   ",if(abs(v-want)<1e-30,"MATCH",if(abs(v+want)<1e-30,"MATCH up to sign","*** DIFFERS ***")));
);
}
quit;
