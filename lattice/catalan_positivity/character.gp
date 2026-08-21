\p 800
G=Catalan;
{klat2(X,Y,V,U,M)=my(K=matkerint([X,V,M,0,0,0; Y,U,0,M,0,0])); mathnf(matrix(2,#K[1,],i,j,K[i,j]));}
{
test(n,k,R,ml)=my(D=lcm(vector(6*n,i,i)),S=D^2,T=2^floor(k*n),M=S*T,
   zr=zudrow(n),nr=nestrow(n),X=zr[1],Y=zr[2],V=nr[1],U=nr[2]);
 my(B=klat2(X,Y,V,U,M), LZ=(X*G-Y)/M, LN=(V*G-U)/M);
 printf("n=%d  sign(Lam_Z)=%d sign(Lam_N)=%d\n",n,sign(LZ),sign(LN));
 for(t=1,#ml, my(l=ml[t], H=Map(), both=0, tot=0, cls=0);
  for(i=-R,R, for(j=-R,R, if(i||j,
    my(cz=i*B[1,1]+j*B[1,2], cn=i*B[2,1]+j*B[2,2], s=sign(cz*LZ+cn*LN));
    my(key=[cz%l, cn%l], cur=if(mapisdefined(H,key), mapget(H,key), [0,0]));
    if(s>0, cur[1]++, if(s<0, cur[2]++)); mapput(H,key,cur))));
  my(ks=Vec(H));
  for(u=1,#ks, my(c=mapget(H,ks[u])); cls++; if(c[1]>0 && c[2]>0, both++));
  printf("   mod %3d : %3d realized classes, %3d of them contain BOTH signs (%.0f%%)\n",
     l, cls, both, 100.0*both/cls));
}
test(12,22.4,60,[3,5,7,16,64]);
test(20,22.4,60,[3,5,7,16,64]);
\q
