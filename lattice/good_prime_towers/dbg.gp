T=towerrun("G",5,[1],4,60);
K0=T[5];recB=T[9];
for(s=0,3, n=5^s; np=5^(s+1); U=125*pratio(mapget(recB,np),mapget(recB,n),5,K0); print(s," U=",U+O(5^6), " v(U-1)=", valuation(U-1,5)));
quit
