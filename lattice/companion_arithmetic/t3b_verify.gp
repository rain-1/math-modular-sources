default(parisize,"6G");
read("lib.gp");
NA = 10000;
NB = 6000;
lg(n,p) = valuation(lcm(1,1),p) + logint(n,p);
{ chk(nm, i, p, fa, fb, NAx, NBx) =
  my(R=ROWS[i], AB=genrow(R,max(NAx,NBx)), A=AB[1], B=AB[2], ba=0, bb=0);
  for(n=1,NAx, if(valuation(A[n+1],p) != fa(n), ba=n; break));
  for(n=1,NBx, if(valuation(B[n+1],p) != fb(n), bb=n; break));
  print(nm, " : a-law holds to n<=", NAx, ": ", if(ba==0,"YES",Str("NO first n=",ba)), "   b-law holds to n<=", NBx, ": ", if(bb==0,"YES",Str("NO first n=",bb)));
}
fa_B(n) = sumdigits(n,3);
fb_B(n) = sumdigits(n,3)-1;
fa_E(n) = 2*sumdigits(n,2);
fb_E(n) = 2*sumdigits(n,2)-2;
fa_eps(n) = 3*sumdigits(n,2) - (n%2);
fb_eps(n) = 3*sumdigits(n,2) - (n%2) - 4;
fa_s(n) = sumdigits(n,2);
fb_s(n) = sumdigits(n,2) - 2*logint(n,2) - 1;
chk("B      p=3, v_3(a_n)=s_3(n), v_3(b_n)=s_3(n)-1        ", 2, 3, fa_B, fb_B, NA, NB);
chk("E      p=2, v_2(a_n)=2s_2(n), v_2(b_n)=2s_2(n)-2      ", 5, 2, fa_E, fb_E, NA, NB);
chk("eps    p=2, v_2(a_n)=3s_2(n)-[n odd], v_2(b_n)=that-4 ", 10, 2, fa_eps, fb_eps, NA, NB);
chk("s10    p=2, v_2(a_n)=s_2(n), v_2(b_n)=s_2(n)-2[lg n]-1", 14, 2, fa_s, fb_s, NA, NB);
chk("s18    p=2, v_2(a_n)=s_2(n), v_2(b_n)=s_2(n)-2[lg n]-1", 15, 2, fa_s, fb_s, NA, NB);
quit;
