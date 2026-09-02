default(parisize,"4G");
read("lib.gp");
read("src.gp");
N = 3000;
{ for(i=1,#ROWS,
  my(R=ROWS[i], r=R[2], AB=genrow(R,N), B=AB[2], bp=badprimes(R), hasC=(i<=12), cf=0);
  if(hasC, cf=CMS[i]);
  for(j=1,#bp,
    my(p=bp[j], neg=List(), nneg=0, mx=0, Rp=0, defect=List(), ndef=0, vv, ee, Rn=0, tail="");
    for(n=1,N,
      vv = valuation(B[n+1],p);
      if(vv<0, nneg++; if(#neg<12, listput(neg,[n,vv])); if(-vv>mx, mx=-vv));
      if(hasC,
        ee = r*valuation(n,p) - if(cf(n)==0, 10^6, valuation(cf(n),p));
        if(ee>Rn, Rn=ee);
        if(max(0,-vv) != Rn, ndef++; if(#defect<8, listput(defect,[n,max(0,-vv),Rn])))));
    if(hasC, tail = Str(" ; v_p(R_N)=", Rn, " ; #n with v_p(den b_n) != v_p(R_n): ", ndef, " ; first (n,v_p(den),v_p(R_n)): ", Vec(defect)));
    print(R[1], " p=", p, " : #{n<=", N, ": v_p(b_n)<0} = ", nneg, " ; max -v_p(b_n) = ", mx, " ; first (n,v_p(b_n)): ", Vec(neg), tail)));
}
quit;
