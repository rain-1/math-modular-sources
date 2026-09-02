default(parisize,"8G");
default(linewrap,0);
BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
read(concat(BASE,"rec_23_10.gp"));
print("(23,10) content of each P_j: ",vector(24,i,content(PV[i])));
read(concat(BASE,"rec_19_13.gp"));
print("(19,13) content of each P_j: ",vector(20,i,content(PV[i])));
print("(19,13) P_19 = (n-18)^7*Q, Q irreducible? ",#factor(PV[20]/(n-18)^7)~);
quit;
