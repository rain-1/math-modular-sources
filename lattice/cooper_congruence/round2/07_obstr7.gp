\\ 07_obstr7.gp -- the obstruction: beta_{s7} is not a combination of twisted CM traces
\\ of modular functions on X_0(7) at discriminant -3m^2.
\\ Basis of weight-0 functions with poles of order <=1 at the two cusps: 1, 1/u, u.
\\ Fricke W_7: u -> 1/(49u).  Anti-invariant: H = 1/u - 49u (trace in i sqrt3 Z).
\\ Invariant:   J = 1/u + 49u (trace real),  and 1 (trace = "class number", size O(m)).
read("lib.gp"); read("heeg.gp");
default(realprecision, 100);
N=7; D0=-3;
{ tr(m) = my(d=-3*m^2, bt, RF, rep, al, ch, om, w, tH=0., tJ=0., t1=0.);
  bt = (5*m)%14;
  RF = redforms(d);
  for(i=1,#RF,
    rep = heegrep(RF[i],N,bt,60);
    if(rep==0, next);
    ch = genchar(rep,D0); om = omeg(rep);
    al = (-rep[2] + I*m*sqrt(3))/(2*rep[1]);
    w = (eta(7*al,1)/eta(al,1))^4;
    tH += ch*(1/w-49*w)/om;
    tJ += ch*(1/w+49*w)/om;
    t1 += ch*1./om);
  [tH/(I*sqrt(3)), tJ, t1];
}
bet = read("beta_s7.txt");
kap = 14/(2*Pi*sqrt(3));
M = 34;
V = vector(M, m, tr(m));
{
print("m  tH*  tJ  t1   beta   beta/tH*   m*beta/((m-kappa)*tH*)");
for(m=1,M,
  if(m%7==0, next);
  print(m,"  ",round(real(V[m][1])),"  ",V[m][2],"  ",V[m][3],"  ",bet[m],"  ",bet[m]/real(V[m][1]),"  ", m*bet[m]/((m-kap)*real(V[m][1])));
);
}
\\ exact-fit test: solve beta = l1*tH* + l2*tJ + l3*t1 on three m, test on the rest
{
my(Mx, rhs, sol, mm=[22,23,24]);
Mx = matrix(3,3,i,j, if(j==1, real(V[mm[i]][1]), if(j==2, real(V[mm[i]][2]), real(V[mm[i]][3]))));
rhs = [bet[mm[1]],bet[mm[2]],bet[mm[3]]]~;
sol = matsolve(Mx,rhs);
print("fit lambda = ", sol~);
print("check (m, beta, predicted):");
for(m=25,M, if(m%7==0,next); print("  ",m,"  ",bet[m],"  ", sol[1]*real(V[m][1])+sol[2]*real(V[m][2])+sol[3]*real(V[m][3])));
}
quit;
