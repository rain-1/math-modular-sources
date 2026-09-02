\\ 22_cprime3.gp -- corrected s18 form of the c'(m) mod 2 law.
read("lib.gp");
CP = [read("20b_cp_s7.txt"), read("20b_cp_s10.txt"), read("20b_cp_s18.txt")];
N  = #CP[1];
print("N = ", N);
{
my(bad);
bad = select(m->my(t=m); while(t%2==0,t=t/2); (CP[3][m]%2==1) != (m%3!=0 && issquare(t)), vector(N,i,i));
print("  s18: c'(m) odd <=> 3 nmid m and odd part of m is a square : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(m->my(t=m); while(t%2==0,t=t/2); while(t%7==0,t=t/7); (CP[1][m]%2==1) != issquare(t), vector(N,i,i));
print("  s7 : c'(m) odd <=> the 14-free part of m is a square      : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
bad = select(m->my(t=m); while(t%2==0,t=t/2); while(t%5==0,t=t/5); (CP[2][m]%2==1) != issquare(t), vector(N,i,i));
print("  s10: c'(m) odd <=> the 10-free part of m is a square      : ", if(#bad==0,"HOLDS",concat("FAILS at ",bad[1])));
}
quit;
