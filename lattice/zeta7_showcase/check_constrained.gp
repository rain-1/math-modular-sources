default(parisize,"16G");
read("/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/lift_lib.gp");
T=(L+1)^7*(L^2-14*L+1)^6;
print("T = ",T);
tv=vector(20,i,polcoeff(T,20-i));
print("t_j j=0..19: ",tv);
mkex(r,D,tv)=matrix(r,(r+1)*(D+1),a,b,my(v=0);if(b==(a+1)*(D+1),v=1);if(b==D+1,v=v-tv[a+1]);v);
EXROWS=mkex(19,13,tv);
print("EXROWS size ",matsize(EXROWS));
p=2^61-1;
setC(p);
res=canon(19,13,p);
print("constrained (19,13) dim = ",matsize(res[1])[1],"  pivots ",res[2]);
{for(ab=0,9, my(a=2*ab+1); if(a<=19, my(b=(19-a)/2, TT=(L+1)^a*(L^2-14*L+1)^b, tw=vector(20,i,polcoeff(TT,20-i))); EXROWS=mkex(19,13,tw); my(rz=canon(19,13,p)); print("  (a,b)=(",a,",",b,") dim=",if(type(rz[1])=="t_MAT",matsize(rz[1])[1],0))));}
quit;
