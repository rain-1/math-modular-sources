default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/01_lib.gp")

PX = 830;   \\ order in x
PY = 412;   \\ order in y  (needs 2*(PY-1) <= PX-1)
AA = 1; BB = -3; CC = 5;   \\ hypothesised rational relation (a,b,c), default of indep_check2.py

N0 = 100; N1 = 400;

\\ ---------- report one function ----------
{
report(nm, cv, pp) =
  my(nn, vv = vector(N1+1, ii, "oo"), t, lo, hi, cnt=0, sx=0,sy=0,sxx=0,sxy=0, sl, wide, v0, v1, prof);
  for(nn = 0, min(N1, #cv-1),
    t = cv[nn+1];
    if(t != 0, vv[nn+1] = valuation(t, pp));
  );
  \\ least squares + min/max of v/n on [N0,N1]
  lo = ""; hi = "";
  for(nn = N0, min(N1, #cv-1),
    t = vv[nn+1];
    if(type(t) == "t_INT",
      cnt++; sx += nn; sy += t; sxx += nn^2; sxy += nn*t;
      if(lo == "" || t/nn < lo, lo = t/nn);
      if(hi == "" || t/nn > hi, hi = t/nn);
    );
  );
  sl = if(cnt >= 2, (cnt*sxy - sx*sy)/(cnt*sxx - sx^2)*1.0, "n/a");
  v0 = vv[N0+1]; v1 = if(N1 < #cv, vv[N1+1], "oo");
  wide = if(type(v0)=="t_INT" && type(v1)=="t_INT", (v1-v0)/300.*1.0, "n/a");
  prof = concat(["["], concat(concat(vector(8, ii, Str(vv[50*ii+1], " ")), ["]"])));
  print("  p=", pp, "  v(c",N0,")=", v0, "  v(c",N1,")=", v1,
        "  wide=", if(type(wide)=="t_REAL", strprintf("%+.5f",wide), wide),
        "  LS=", if(type(sl)=="t_REAL", strprintf("%+.5f",sl), sl),
        "  min v/n=", if(lo=="", "-", strprintf("%+.5f", lo*1.0)),
        "  max v/n=", if(hi=="", "-", strprintf("%+.5f", hi*1.0)),
        "  nonzero=", cnt, "/", N1-N0+1);
  print("        v_",pp," at n=50,100,...,400: ", concat(vector(8, ii, Str(vv[50*ii+1], " "))));
}

\\ ---------- build ----------
gettime();
HA = solveode(0,0,PX); HB = solveode(1,0,PX); HC = solveode(0,1,PX);
if(HA[1]!=1 || HA[2]!=3 || HA[3]!=15 || HA[4]!=93, error("HA check failed"));
if(HB[3] != 23/4, error("HB check failed"));
if(HC[4] != 343/9, error("HC check failed"));

H  = vector(PX, ii, AA*HA[ii] + BB*HB[ii] + CC*HC[ii]);
DD = lcm(vector(PX, ii, denominator(H[ii])));
Hi = vector(PX, ii, H[ii]*DD);
Hw = composew(Hi, PX);
Gi = toy(vector(PX, ii, Hi[ii]+Hw[ii]), PX, PY);
G  = vector(PY, ii, Gi[ii]/DD);
if(G[1] != 2 || G[2] != -101/2 || G[3] != 15285/8, error("G check failed"));

\\ B4 = to_y( anti * (Li2(x) - Li2(w(x))) ),  anti = x - w(x)
D2 = lcm(vector(PX, ii, ii))^2;
Li2 = vector(PX); for(nn=1, PX-1, Li2[nn+1] = D2/nn^2);
Li2w = composew(Li2, PX);
B4i = toy(antimul(vector(PX, ii, Li2[ii]-Li2w[ii]), PX), PX, PY);
B4 = vector(PY, ii, B4i[ii]/D2);
if(B4[2] != -4 || B4[3] != 4/9 || B4[4] != 31/900, error("B4 check failed"));

B1 = vector(PY); B1[1] = 1;
B2 = vector(PY); for(nn=2, PY-1, B2[nn+1] = 2*(nn-2)!*nn!/(2*nn)!);
B3 = vector(PY); for(nn=1, PY-1, B3[nn+1] = ((nn-1)!)^2/(2*nn)!);
B5 = vector(PY); for(nn=1, PY-1, B5[nn+1] = ((nn-1)!)^2/((2*nn-1)!*(2*nn-1)));
B6 = vector(PY); for(nn=1, PY-1, B6[nn+1] = ((nn-1)!)^2/(nn*(2*nn)!));

\\ literal ports of indep_check2.py lines 91,96,97
B7   = iy(concat([0], vector(PY-1, nn, B4[nn+1]/nn)));
IG1  = iy(concat([0], vector(PY-1, nn, G[nn+1]/nn)));
IG2  = iy(concat([0,0], vector(PY-2, nn, G[nn+2]/nn)));
\\ "intended" (un-shifted) variants, for comparison
B7b  = concat([0], vector(PY-1, nn, B4[nn+1]/nn^2));
IG1b = concat([0], vector(PY-1, nn, G[nn+1]/nn^2));
IG2b = concat([0,0], vector(PY-2, nn, G[nn+2]/((nn+1)*nn^2)));

{
fns = [["B1",B1],["B2",B2],["B3",B3],["B4",B4],["B5",B5],
       ["G",G],["G'",dy(G)],["G''",dy(dy(G))],["G'''",dy(dy(dy(G)))],
       ["B6",B6],["B7",B7],["intG",iy(G)],
       ["int(G-G0)/y",IG1],["int(G-G0-G1y)/y^2",IG2],
       ["[alt]B7=B4n/n^2",B7b],["[alt]int(G-G0)/y",IG1b],["[alt]int(G-G0-G1y)/y^2",IG2b]];
}


print("build time ", gettime(), " ms;  PX=",PX," PY=",PY," (a,b,c)=(",AA,",",BB,",",CC,")");
print("digits(DD)=", #digits(DD), "  v_2(DD)=", valuation(DD,2), "  v_3(DD)=", valuation(DD,3));
print("");
print("=== PART 1: y-coefficients of the 14 CDT functions (window n in [",N0,",",N1,"]) ===");
{
foreach(fns, fn,
  print(fn[1], "   (length ", #fn[2], ")");
  report(fn[1], fn[2], 2);
  report(fn[1], fn[2], 3);
  print("");
);
}

print("=== PART 2: x-coefficients of H_A, H_B, H_C ===");
{
foreach([["H_A",HA],["H_B",HB],["H_C",HC],["H=aA+bB+cC",H]], fn,
  print(fn[1]);
  report(fn[1], fn[2], 2);
  report(fn[1], fn[2], 3);
  print("");
);
}
