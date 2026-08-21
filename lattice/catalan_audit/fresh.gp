\p 800
default(parisize,4000000000);
G=Catalan;
rd(f)=my(v=readstr(f),h=Map()); for(i=1,#v, my(a=Vec(Str(v[i]))); ); v;
Zl=readstr("zud_exact.txt"); Nl=readstr("nest_exact.txt");
{
getrow(lines,n)=my(r); for(i=1,#lines, my(w=Vec(strsplit(lines[i]," "))); if(#w==3 && eval(w[1])==n, r=[eval(w[2]),eval(w[3])])); r;
}
{
for(n=6,13,
 my(XY=getrow(Zl,n), VU=getrow(Nl,n), X=XY[1],Y=XY[2],V=VU[1],U=VU[2]);
 my(S=lcm(vector(6*n,i,i))^2);
 if(X%S||V%S, print("S nodiv"); next);
 my(a1=X/S, a2=V/S, h=a1*U-a2*Y, vh=valuation(h,2), T=2^vh, M=S*T);
 printf("n=%d v2(h)=%d (%.3f bits/n) log2 S/n=%.3f sigma=%.3f\n",n,vh,vh/n,log(S)/log(2)/n,log(M)/n);
 /* HNF basis of K = {c : a1 c1 + a2 c2 = 0 mod T, Y c1 + U c2 = 0 mod M} */
 my(Mt=matconcat([[a1,a2;Y,U],matdiagonal([T,M])]));
 /* kernel mod: use matrix over Z of the map Z^2 -> Z/T x Z/M ; lattice K = preimage */
 my(L=matrix(2,4)); 
 /* build generators of K as columns: solve via HNF of the lattice {(c1,c2,k1,k2): a1c1+a2c2=k1 T, Yc1+Uc2=k2 M} */
 my(B=matkerint(matconcat([[a1,a2,-T,0;Y,U,0,-M]])));
 my(gens=vector(#B[1,],j,[B[1,j],B[2,j]]));
 my(K=mathnf(Mat(vector(#gens,j,[gens[j][1],gens[j][2]]~))));
 my(b1=[K[1,1],K[2,1]], b2=[K[1,2],K[2,2]], covol=abs(matdet(K)));
 printf("   covolK=%s  M=%s  covol==M? %d\n", covol, M, covol==M);
 my(E1=log(abs(X*G-Y))/n, E2=log(abs(V*G-U))/n, sig=log(M)/n, x=(sig+E2-E1)/2);
 my(w=floor(exp((2*x-sig)*n))); if(w<1,w=1);
 
 my(Kw=[b1[1],b2[1];b1[2]*w,b2[2]*w], red=Kw*qflll(Kw));
 for(j=1,2, my(c1=red[1,j], c2w=red[2,j]);
   if(c2w%w!=0, printf("   cand %d: not exact lattice pt\n",j); next);
   my(c2=c2w/w, nq=c1*X+c2*V, np=c1*Y+c2*U);
   if(nq%M||np%M, printf("   cand %d: NOT divisible by M\n",j); next);
   my(q=nq/M, p=np/M);
   if(q==0, printf("   cand %d: q=0\n",j); next);
   printf("   cand %d: log|q|/n=%.4f  log|qG-p|/n=%.4f  gcd=%d\n", j, log(abs(q))/n, log(abs(q*G-p))/n, gcd(q,p)));
);
}
\q
