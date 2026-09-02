\\ 25_mod9.gp -- the mod-9 pair law for s18 suggested by the period-18 pattern.
read("lib.gp");
G  = [read("20c_gamma_s7.txt"), read("20c_gamma_s10.txt"), read("20c_gamma_s18.txt")];
N  = #G[1];
print("N = ", N);
{
my(bad);
bad = select(n->(n%3==1 && (G[3][n]*G[3][n+1]+1)%9!=0), vector(N-1,i,i));
print("  s18: gamma(n)gamma(n+1) = -1 mod 9 for n = 1 mod 3 : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(n->(n%3==1 && (G[3][n]*G[3][n+1]+1)%27!=0), vector(N-1,i,i));
print("  ... mod 27 : ", if(#bad==0,"HOLDS",concat("FAILS at n=",bad[1])));
bad = select(n->((G[3][n]-G[3][if(n<=18,n+18,18-n+18*((n-1)\18+1))])%9!=0 && n%3!=0), vector(18,i,i));
print("  s18: gamma(n) = gamma(18-n) mod 9 for 3 nmid n, n<18 : ",
      if(#select(n->(n%3!=0 && n<18 && (G[3][n]-G[3][18-n])%9!=0), vector(17,i,i))==0,"HOLDS","FAILS"));
bad = select(n->(n%3!=0 && (G[1][n]*G[1][n+1])%1!=0), vector(N-1,i,i));
}
\\ analogue for the other rows / other moduli?
{
my(bad);
bad = select(n->(n%2==1 && n%7!=0 && (n+1)%7!=0 && (G[1][n]*G[1][n+1])%4!=0), vector(N-1,i,i));
print("  s7 : gamma(n)gamma(n+1) = 0 mod 4 for n odd, 7 nmid n(n+1) : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
}
quit;
