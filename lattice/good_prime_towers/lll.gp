/* lll.gp -- identification tests on the tower limits Lambda_a.
   Input: LAM lines (row p a w chi vLam ndig U) collected from the run logs.
   Reports best relation height against the heuristic noise floor p^(ndig/m). */

AP4 = Map();  /* a_p of eta(2z)^4 eta(4z)^4 in S_4(Gamma_0(8)) -- Apery zeta(3) */
mapput(AP4,5,-2); mapput(AP4,7,24); mapput(AP4,11,-44); mapput(AP4,13,22);
AP3 = Map();  /* a_p of eta(2z)^3 eta(6z)^3 in S_3(Gamma_0(12),chi_-3) */
mapput(AP3,5,0); mapput(AP3,7,2); mapput(AP3,11,0); mapput(AP3,13,-22);

/* unit root of X^2 - ap X + p^k, if ordinary; 0 if not */
{
unitroot(ap,p,k,pr) =
  my(d,r);
  if(ap%p==0, return(0));
  r = ap + O(p^pr);
  /* Newton on x^2-ap x+p^k, starting from ap */
  for(i=1,pr+5, r = r - (r^2-ap*r+p^k)/(2*r-ap));
  r;
}

{
hgt(v) = vecmax(apply(abs,Vec(v)));
}

{
test(lab, vec, K, p) =
  my(r,h,fl);
  r = lindep(vec);
  h = hgt(r);
  fl = p^(K/#vec);
  print("    ",lab,": height ",h,"   floor ",floor(fl),
        if(h < fl/4, "   <-- BELOW FLOOR", ""));
}
