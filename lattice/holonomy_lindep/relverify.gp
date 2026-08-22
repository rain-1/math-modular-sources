default(parisizemax, 8G);
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/ops.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/sc_rows.gp");
mkcf(op) = { my(dz=op[3], Ps=op[4]); vector(dz+1, i, subst(Ps[i], 'X, 'n-(i-1))); };
findop(cl) = { for(t=1,#OPS, if(OPS[t][2]==cl, return(OPS[t]))); 0; };
prim(v) = { my(g=0); for(i=1,#v, g=gcd(g,v[i])); if(g==0, return(v)); v=v/g;
            my(i=1); while(i<=#v && v[i]==0, i++); if(i<=#v && v[i]<0, v=-v); v; };
lidep(vals, D) = {
  my(dl = max(20, floor(0.45*D)), dh = max(30, floor(0.85*D)), pv, ph, sv);
  sv = default(realprecision);
  default(realprecision, dl); pv = prim(Vec(lindep(vector(#vals,i,vals[i]*1.0))));
  default(realprecision, dh); ph = prim(Vec(lindep(vector(#vals,i,vals[i]*1.0))));
  default(realprecision, sv);
  if(pv == ph, ph, 0);
};
N=800; NB=750;
\\ full relation lattice of [1,xi_1..xi_{dz-1}] by iterated elimination
{
foreach(["4.8.1","4.5.101","4.8.2","4.8.6"], cl,
  my(op=findop(cl), dz=op[3], cf=mkcf(op));
  my(a=genseq(cf,dz,[1],N));
  my(xs=vector(dz-1,j,compan(cf,dz,j,N)));
  default(realprecision,1600);
  my(xi=vector(dz-1), cert=vector(dz-1));
  for(j=1,dz-1,
    my(rN=(xs[j][N+1]*1.0)/a[N+1], rB=(xs[j][NB+1]*1.0)/a[NB+1]);
    xi[j]=rN; my(d=abs(rN-rB)); cert[j]=if(d==0,1500,min(1500,floor(-log(d)/log(10)))));
  my(D=vecmin(cert));
  print("### ", op[1], " ", cl, " D=", D);
  my(vals = concat([1], vector(dz-1,j,xi[j])));
  my(labs = concat(["1"], vector(dz-1,j,Str("xi",j))));
  my(nrel = 0);
  while(1,
    my(r = lidep(vals, D));
    if(r === 0, break);
    default(realprecision, min(D,400));
    my(res = abs(sum(i=1,#vals, r[i]*vals[i])));
    my(H = vecmax(vector(#r,i,abs(r[i]))));
    print("  RELATION #", nrel+1, ": ", r, " on ", labs, "  height=", H,
          "  log10|res|=", if(res==0,-999,round(log(res)/log(10))), "  (working D=",min(D,400),")");
    default(realprecision,1600);
    \\ eliminate: drop a coordinate with nonzero coeff
    my(p=0); for(i=1,#r, if(r[i]!=0, p=i; break));
    my(keep = setminus(Set(vector(#vals,i,i)), Set([p])));
    vals = vector(#keep,i,vals[keep[i]]); labs = vector(#keep,i,labs[keep[i]]);
    nrel++;
    if(#vals <= 1, break);
  );
  print("  #independent relations = ", nrel, "  => prk = ", (dz-1) - nrel);
  print("");
);
}
quit;
