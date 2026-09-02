default(parisize,1000000000);
c=read("cn400.txt"); d=read("dn.txt");
cp=vector(8,i,c[i]+if(i>1,c[i-1],0)); dp=vector(8,i,d[i]+if(i>1,d[i-1],0));
print("c'_0..7 = ",cp);
print("d'_1..4 = ",vector(4,i,dp[i+1]));
print("(1728/209) d'_n/c'_n, n=1..4: ",vector(4,i,1728/209*dp[i+1]/cp[i+1]));
print("K7' = ",(8-4*sqrt(3))*3^(27/4)*(739-356*sqrt(3))/2^(77/6)*gamma(1/3)^12/Pi^(19/2));
print("7alpha'/K7' = ",(14161/720)/((8-4*sqrt(3))*3^(27/4)*(739-356*sqrt(3))/2^(77/6)*gamma(1/3)^12/Pi^(19/2)));
quit;
