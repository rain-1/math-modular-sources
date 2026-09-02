\\ 24_lucas.gp -- Lucas / Dwork / periodicity tests for gamma mod p.
read("lib.gp");
G  = [read("20_gamma_s7.txt"), read("20_gamma_s10.txt"), read("20_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
PL = vector(14,i,prime(i));

print();
print("=== (a) gamma(pn) = gamma(p)gamma(n) mod p   [fail/total] ===");
{
for(k=1,3, my(s="");
  for(ip=1,14, my(p=PL[ip], f=0, t=0, gp1=G[k][p]);
    for(n=1,N\p, t++; if((G[k][p*n]-gp1*G[k][n])%p!=0, f++));
    s=concat(s,concat(concat(concat(concat("  p=",p),": "),concat(f,"/")),t)));
  print(NAM[k], s));
}
print();
print("=== (c) gamma(pn) = psi(p)gamma(n) mod p   [fail/total] ===");
{
for(k=1,3, my(s="");
  for(ip=1,14, my(p=PL[ip], f=0, t=0, ps=psival(k,p));
    for(n=1,N\p, t++; if((G[k][p*n]-ps*G[k][n])%p!=0, f++));
    s=concat(s,concat(concat(concat(concat("  p=",p),": "),concat(f,"/")),t)));
  print(NAM[k], s));
}
print();
print("=== (c') best constant: is gamma(pn)/gamma(n) mod p constant in n? ===");
{
for(k=1,3,
  for(ip=1,14, my(p=PL[ip], rs=List(), t=0);
    for(n=1,N\p, if(G[k][n]%p!=0, t++; listput(rs, lift(Mod(G[k][p*n],p)/Mod(G[k][n],p)))));
    rs=Set(Vec(rs));
    print("  ",NAM[k]," p=",p,"  #distinct ratios = ",#rs,"  over ",t," usable n",
          if(#rs<=3, concat("   values: ",rs), "")));
);
}
print();
print("=== (c'') mod p^2 : gamma(pn) = psi(p) gamma(n) mod p^2 [fail/total] ===");
{
for(k=1,3, my(s="");
  for(ip=1,14, my(p=PL[ip], f=0, t=0, ps=psival(k,p));
    for(n=1,N\p, t++; if((G[k][p*n]-ps*G[k][n])%(p^2)!=0, f++));
    s=concat(s,concat(concat(concat(concat("  p=",p),": "),concat(f,"/")),t)));
  print(NAM[k], s));
}
print();
print("=== (b) gamma(pn+r) = gamma(n)gamma(r) mod p, 1<=r<p, n>=1 [fail/total] ===");
{
for(k=1,3, my(s="");
  for(ip=1,14, my(p=PL[ip], f=0, t=0);
    for(n=1,(N-1)\p, for(r=1,p-1, if(p*n+r<=N, t++;
      if((G[k][p*n+r]-G[k][n]*G[k][r])%p!=0, f++))));
    s=concat(s,concat(concat(concat(concat("  p=",p),": "),concat(f,"/")),t)));
  print(NAM[k], s));
}
print();
print("=== (d) eventual periodicity of gamma(n) mod p: smallest (offset<=40, period<=400) ===");
{
for(k=1,3,
  for(ip=1,14, my(p=PL[ip], found=0);
    for(off=0,40, if(found,break);
      for(T=1,400, if(off+2*T>N, break);
        my(ok=1);
        for(n=off+1,N-T, if((G[k][n]-G[k][n+T])%p!=0, ok=0; break));
        if(ok, print("  ",NAM[k]," p=",p,"  PERIODIC offset=",off," period=",T); found=1; break)));
    if(!found, print("  ",NAM[k]," p=",p,"  not periodic (offset<=40, period<=400)")));
);
}
quit;
