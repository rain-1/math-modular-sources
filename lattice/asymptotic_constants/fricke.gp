default(parisize,"10G");
default(realprecision,120);
MQ = 700;
Dop(S) = 'q*deriv(S,'q);
etaq(DL,EL)={
  my(r=1+O('q^(MQ+1)), s=0);
  for(i=1,#DL, s += DL[i]*EL[i]; for(m=1,MQ\DL[i], r *= (1-'q^(DL[i]*m)+O('q^(MQ+1)))^EL[i]));
  'q^(s/24)*r;
}
ev(S,z)=subst(truncate(S),'q,z);
ROWS = List();
listput(ROWS, ["alpha", 12, 10, 9, [1,2,3,4,6,12], [-4,4,4,-4,-4,4]]);
listput(ROWS, ["gamma",  6, 17, 72, [1,2,3,6], [-5,1,-1,5]]);
listput(ROWS, ["eps",    8, 12, 32, [1,2,4,8], [-4,2,-2,4]]);
listput(ROWS, ["zeta",   9, 9, 27, [1,3,9], [-3,0,3]]);
listput(ROWS, ["s7",     7, 13, 49, [1,7], [-4,4]]);
listput(ROWS, ["s10",   10, 6, 25, [1,2,5,10], [-2,-2,2,2]]);
listput(ROWS, ["s18",   18, 14, 1, [1,2,3,6,9,18], [-6,6,12,-12,-6,6]]);
th=1.9;
{
for(i=1,#ROWS,
  my(R=ROWS[i], nm=R[1], N=R[2], BB=R[3], CC=R[4]);
  my(uq=etaq(R[5],R[6]), Fq=Dop(uq)/uq, tq=uq/(1+BB*uq+CC*uq^2));
  my(tau=(cos(th)+I*sin(th))/sqrt(N), taup=-1/(N*tau));
  my(q1=exp(2*Pi*I*tau), q2=exp(2*Pi*I*taup));
  print(nm,": u(W tau) u(tau) = ", ev(uq,q2)*ev(uq,q1), "   1/C = ", 1./CC);
  print("       F|_2 W_N / F  = ", ev(Fq,q2)/((sqrt(N)*tau)^2*ev(Fq,q1)));
  print("       t(W tau)-t(tau) = ", ev(tq,q2)-ev(tq,q1));
);
}
quit;
