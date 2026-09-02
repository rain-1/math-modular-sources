default(parisize,"16G");
BASE="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/zeta7/";
print("18951045120 = ",factor(18951045120));
print("109182976 = ",factor(109182976));
print("207884386304 = ",factor(207884386304));
D4=read(concat(BASE,"dn.txt"));
D9=read(concat(BASE,"dn900.txt"));
print("#dn.txt=",#D4,"  #dn900.txt=",#D9);
ok=1;
{for(i=1,min(#D4,#D9), if(D4[i]!=D9[i],ok=0;break));}
print("dn900 agrees with dn on the first ",min(#D4,#D9)," entries: ",ok);
quit;
