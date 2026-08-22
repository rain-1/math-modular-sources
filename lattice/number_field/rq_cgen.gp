default(parisizemax, 6000000000);
default(realprecision, 120);
Ds = [2,3,5,6,7,10,11,13,14,15,17,19,21,22,23,26,29,30,31,33,34,35,37,38,39,41,42,43,46,47,51,53,55,57,58,59,61,62,65,66,67,69,70,71,73,74,77,78,79,82,83,85,86,87,89,91,93,94,95,97];
{
for(i=1,#Ds,
  my(D=Ds[i], t, s, disc, eps, L, L2, w1, w2, ev);
  if(D%4==1, t=1; s=(D-1)/4; disc=D, t=0; s=D; disc=4*D);
  eps = quadunit(disc);
  w1 = (t+sqrt(t^2+4*s))/2; w2 = (t-sqrt(t^2+4*s))/2;
  L = List();
  for(m=-6,6, my(p=eps^m, p0=real(p), p1=imag(p));
     listput(L,[p0,p1]); listput(L,[-p0,-p1]));
  L2 = List();
  for(c0=-20,20, for(c1=-20,20,
     my(nn=c0^2+t*c0*c1-s*c1^2);
     if(nn!=0 && abs(nn)<=16, listput(L2,[c0,c1]))));
  \\ validation list: rational C with |C|<=130
  my(L3=List()); for(c0=-130,130, if(c0!=0, listput(L3,[c0,0])));
  print("D ",D," ",t," ",s);
  print("NC1 ",#L);
  for(j=1,#L, my(c0=L[j][1],c1=L[j][2]);
     print("C ",c0," ",c1," ",precision(c0+c1*w1,30)," ",precision(c0+c1*w2,30)));
  print("NC2 ",#L2);
  for(j=1,#L2, my(c0=L2[j][1],c1=L2[j][2]);
     print("S ",c0," ",c1," ",precision(c0+c1*w1,30)," ",precision(c0+c1*w2,30)));
  print("NC3 ",#L3);
  for(j=1,#L3, print("Q ",L3[j][1]," 0 ",L3[j][1],".0 ",L3[j][1],".0"));
  print("ENDD");
);
}
quit;
