/* 03_scan.gp -- for every genus-0 host Gamma_0(N)+W with 4|N, enumerate the
   joint Atkin-Lehner eigenspaces of M_1(N,chi_{-4}) and M_3(N,chi_{-4}),
   report their dimensions and the number of free parameters in the pair
   (F, Phi_0)  with  a_0(F)=1,  a_0(Phi_0)=0, a_1(Phi_0)=1.
   Matching eigenvalue is required:  (Phi_0/F)|_2 W_Q = (Phi_0|_3 W_Q)/(F|_1 W_Q).
*/
default(parisizemax, 8000000000);
default(realprecision, 40);

{HOSTS = [ [4,[]], [4,[4]], [8,[]], [8,[8]], [12,[]], [12,[3]], [12,[4]], [12,[3,4]],
          [16,[]], [16,[16]], [20,[4]], [20,[4,5]], [24,[8]], [24,[3,8]],
          [28,[7]], [28,[4,7]], [32,[32]], [36,[4]], [36,[4,9]],
          [44,[4,11]], [56,[7,8]], [60,[3,4,5]], [92,[4,23]] ];}

/* eigenvalues (complex) of W_Q on a space: returns the diagonalised data
   [evals, eigenvector-matrix] assuming W_Q is diagonalisable with evals in
   {+-1,+-i}.  We instead return, for a requested eigenvalue ev, the kernel. */
{eigspace(mfx, qs, evs) =
  my(d = mfdim(mfx));
  if(d==0, return(matrix(0,0)));
  my(pr = matid(d));
  for(j=1,#qs,
    my(al = mfatkininit(mfx,qs[j]), mm = al[2]/al[3]);
    /* projector onto eigenvalue evs[j] : prod over other eigenvalues */
    my(others = select(z->abs(z-evs[j])>1e-10, [1,-1,I,-I]));
    my(p1 = matid(d));
    for(k=1,#others, p1 = p1*(mm - others[k]*matid(d))/(evs[j]-others[k]));
    pr = pr*p1;
  );
  my(rr = matrix(d,d,i,j, bestappr(real(pr[i,j]),10^18) + I*bestappr(imag(pr[i,j]),10^18)));
  matimage(rr);
}

print("N W | ev | dim M1eig dim M3eig | params");
{for(h=1,#HOSTS,
  my(nn=HOSTS[h][1], qs=HOSTS[h][2]);
  my(mf1 = mfinit([nn,1,-4],4), mf3 = mfinit([nn,3,-4],4));
  my(d1 = mfdim(mf1), d3 = mfdim(mf3));
  my(nev = 2^#qs);
  print("--- N=",nn," W=",qs,"   dim M1=",d1,"  dim M3=",d3);
  /* possible eigenvalue for each Q: determined by W_Q^2 = -1 or +1 */
  my(cand = vector(#qs, j, if(qs[j]%4==0 || qs[j]%4==3, [I,-I], [1,-1])));
  my(idx = vector(#qs,j,1));
  for(cc=0, nev-1,
    my(evs = vector(#qs, j, cand[j][ 1 + bittest(cc,j-1) ]));
    my(e1 = eigspace(mf1,qs,evs), e3 = eigspace(mf3,qs,evs));
    my(f = matsize(e1)[2], e = matsize(e3)[2]);
    print("   ev=",evs,"  dimF=",f,"  dimPhi=",e,"  params=", max(0,f-1)+max(0,e-2));
  );
);}
quit;
