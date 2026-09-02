default(realprecision,60);
G=znstar(5,1);
chars=[[1],[2],[3],[0]];
for(j=1,4, my(ch=[j%4]); print("char idx ",j," : ", vector(5,a,chareval(G,ch,a))));
quit;
