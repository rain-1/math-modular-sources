default(parisize,"3G");
read("lib.gp");
read("src.gp");
N = 3000;
{ for(i=1,12,
  my(R=ROWS[i], r=R[2], AB=genrow(R,N), B=AB[2], dn=1, Rn=1, cf=CMS[i]);
  my(mxe=0, small=List(), cnt=0, tot=0, big=0, mid=0, sml=0, fa, T, dv);
  for(n=1,N,
    dn = lcm(dn,n);
    my(t=n^r, g=gcd(t, cf(n))); Rn = lcm(Rn, t/g);
    dv = denominator(B[n+1]);
    T = Rn/dv;
    if(T>1, cnt++; fa=factor(T);
      for(j=1,#fa[,1],
        tot++;
        if(fa[j,2]>mxe, mxe=fa[j,2]);
        if(fa[j,1]<=50, listput(small, [n,fa[j,1],fa[j,2]]); sml++, if(2*fa[j,1]>n, big++, mid++)))));
  print(R[1], " r=", r, ": #{n<=",N,": den<R_n} = ", cnt, " ; total (p,e) defects = ", tot, " ; max defect exponent = ", mxe, " ; defects with p<=50: ", sml, " ; with 50<p<=n/2: ", mid, " ; with p>n/2: ", big);
  if(sml>0, print("     small-prime defects (n,p,e), first 12: ", vector(min(12,#small), j, small[j]))));
}
quit;
