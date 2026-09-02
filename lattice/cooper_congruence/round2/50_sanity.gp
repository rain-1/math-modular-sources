\\ 50_sanity.gp -- fast lib agrees with the old lib and with the stored data files.
read("50_lib.gp");
M = 300;
old = [read("20_cp_s7.txt"), read("20_cp_s10.txt"), read("20_cp_s18.txt")];
oldb= [read("20_beta_s7.txt"), read("20_beta_s10.txt"), read("20_beta_s18.txt")];
gettime();
{
for(k=1,3,
  my(cp=CPvec(k,M), b=Bvec(k,cp), ok=1, okb=1);
  for(m=1,M, if(cp[m]!=old[k][m], ok=0; print("  cp mismatch at m=",m)));
  for(m=1,M, if(b[m]!=oldb[k][m], okb=0));
  print("row ",NAM[k],": c' matches stored: ",ok,"   beta matches stored: ",okb);
);
}
print("time ", gettime(), " ms");
quit;
