/* 05_rational.gp -- which genus-0 AL quotient hosts carry a RATIONAL
   chi_{-4} weight-1 eigenvector (needed for A,B in Q[[t]]).
   W_Q acts on M_k(N,chi_{-4}) with W_Q^2 = -1 (4|Q or Q=3 mod 4) or +1
   (Q = 1 mod 4).  The eigenprojector is rational iff the eigenspace is. */
default(parisizemax, 8000000000);
default(realprecision, 60);

{HOSTS = [ [4,[4]], [8,[8]], [12,[3]], [12,[4]], [12,[3,4]], [16,[16]],
          [20,[4]], [20,[4,5]], [24,[8]], [24,[3,8]], [28,[7]], [28,[4,7]],
          [32,[32]], [36,[4]], [36,[4,9]], [44,[4,11]], [56,[7,8]],
          [60,[3,4,5]], [92,[4,23]] ];}

{ratproj(mm, ev) =
  my(d=matsize(mm)[1]);
  my(others = select(z->abs(z-ev)>1e-10, [1,-1,I,-I]));
  my(p1 = matid(d));
  for(k=1,#others, p1 = p1*(mm - others[k]*matid(d))/(ev-others[k]));
  my(ok=1, rr=matrix(d,d));
  for(i=1,d, for(j=1,d,
    my(z=p1[i,j]);
    if(abs(imag(z))>1e-25, ok=0);
    my(a=bestappr(real(z),10^6));
    if(abs(a-real(z))>1e-25, ok=0);
    rr[i,j]=a;
  ));
  [ok, rr];
}

print("N  W        ev        M1: rational? dim   M3: rational? dim");
{for(h=1,#HOSTS,
  my(nn=HOSTS[h][1], qs=HOSTS[h][2]);
  my(mf1=mfinit([nn,1,-4],4), mf3=mfinit([nn,3,-4],4));
  my(d1=mfdim(mf1), d3=mfdim(mf3));
  my(m1s = vector(#qs,j, my(al=mfatkininit(mf1,qs[j])); al[2]/al[3]));
  my(m3s = vector(#qs,j, my(al=mfatkininit(mf3,qs[j])); al[2]/al[3]));
  my(cand = vector(#qs, j, if(qs[j]%4==1, [1,-1], [I,-I])));
  for(cc=0, 2^#qs-1,
    my(evs = vector(#qs, j, cand[j][ 1 + bittest(cc,j-1) ]));
    my(ok1=1, p1=matid(d1), ok3=1, p3=matid(d3));
    for(j=1,#qs,
      my(r=ratproj(m1s[j],evs[j])); ok1 = ok1 && r[1]; if(r[1], p1 = p1*r[2]);
      my(s=ratproj(m3s[j],evs[j])); ok3 = ok3 && s[1]; if(s[1], p3 = p3*s[2]);
    );
    my(dd1 = if(ok1, matrank(p1), -1), dd3 = if(ok3, matrank(p3), -1));
    print(nn," ",qs,"\t",evs,"\t", if(ok1,Str("Q dim=",dd1),"irrational"), "\t",
          if(ok3,Str("Q dim=",dd3),"irrational"));
  );
);}
quit;
