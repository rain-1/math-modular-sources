read("maass.gp");
default(realprecision,50);
FS = vector(3,k,mkfser(k,600));
{
for(k=1,3,
  my(N=LEV[k], t=0.113+0.317*I, A, B);
  A = fhatC(k,t); B = fhatQ(t,N,FS[k],550);
  print("row ",NAM[k],"  fhat closed = ",A[2]);
  print("            fhat q-ser  = ",B);
);
}
quit;
