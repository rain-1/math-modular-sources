default(parisize,"20G");
read("/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/analyze_lib.gp");
read(concat(BASE,"rec_19_13_raw.gp"));
RR=19;
DD=13;
v=prim(BVL[1]);
PV=tovec2pol(v,RR,DD);
print("== shape (19,13) canonical ==");
print("max coeff digits over all P_j: ",maxdig(v));
print("digits of P_j leading coeff and content:");
{for(jj=0,RR, print("  j=",jj,"  deg=",poldegree(PV[jj+1],n),"  content digits=",#Str(numerator(content(PV[jj+1]))),"  maxcoefdig=",maxdig(Vec(PV[jj+1]))));}
ch=chiof(v,RR,DD);
print("chi content = ",content(ch));
print("chi/content factored: ");
print(factor(ch/content(ch)));
print("chi equals c*(L+1)^7*(L^2-14L+1)^6 ? ",ch/content(ch)==(L+1)^7*(L^2-14*L+1)^6);
print("");
print("symmetry test:");
{for(A=-40,60, my(w=reflect(v,RR,DD,A)); if(w==v,print("  A=",A," eps=+1")); if(w==-v,print("  A=",A," eps=-1")));}
print("symmetry scan done");
quit;
