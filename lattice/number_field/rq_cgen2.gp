default(parisizemax, 6000000000);
default(realprecision, 120);
Ds = [2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,39,41,42,43,46,47,51,53,55,57,58,59,61,62,65,66,67,69,70,71,73,74,77,78,79,82,83,85,86,87,89,91,93,94,95,97];
{
for(i=1,#Ds,
  my(D=Ds[i], t, s, disc, eps, w1, w2, U, T, seen);
  if(D%4==1, t=1; s=(D-1)/4; disc=D, t=0; s=D; disc=4*D);
  eps = quadunit(disc);
  w1=(t+sqrt(t^2+4*s))/2; w2=(t-sqrt(t^2+4*s))/2;
  \\ U: units modulo squares of units: {1,-1,eps,-eps}
  U = [[1,0],[-1,0],[real(eps),imag(eps)],[-real(eps),-imag(eps)]];
  \\ T: C with 1<=|N(C)|<=16, |c0|,|c1|<=20, deduped modulo multiplication by eps^2
  T = List(); seen = Map();
  my(ww=Mod(y,y^2-t*y-s), ee=real(eps)+imag(eps)*ww);
  for(c0=-20,20, for(c1=-20,20,
    my(nn=c0^2+t*c0*c1-s*c1^2);
    if(nn!=0 && abs(nn)<=16,
      my(z=c0+c1*ww, key=0, found=0);
      for(m=-8,8, my(zz=lift(z*ee^(2*m)), k2=[polcoeff(zz,0),polcoeff(zz,1)]);
        if(mapisdefined(seen,k2), found=1; break));
      if(!found, mapput(seen,[c0,c1],1); listput(T,[c0,c1])))));
  print("D ",D," ",t," ",s);
  print("NU ",#U);
  for(j=1,#U, my(c0=U[j][1],c1=U[j][2]); print("U ",c0," ",c1," ",precision(c0+c1*w1,30)," ",precision(c0+c1*w2,30)));
  print("NT ",#T);
  for(j=1,#T, my(c0=T[j][1],c1=T[j][2]); print("T ",c0," ",c1," ",precision(c0+c1*w1,30)," ",precision(c0+c1*w2,30)));
  print("ENDD");
);
}
quit;
