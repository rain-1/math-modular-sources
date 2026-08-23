/* 07_ident.gp -- identification battery for the Apery limits xi of four-term rows.
 * usage:  gp -q -s 2000000000 07_ident.gp     reading  out/xi.txt  (lines "label xi")
 * Avoids the PARI builtins psi, M, Phi, S, cmp.
 */
default(realprecision, 90);
PREC = 60;

NAMES = List(); VALS = List();
addc(nm, v) = { listput(NAMES, nm); listput(VALS, v); }

{
  addc("zeta(2)", zeta(2));
  foreach([-3,-4,-7,-8,-11,-15,-19,-20,-23,-24,5,8,12,13,17,21,24], dd,
          addc(Str("L(2,chi_", dd, ")"), lfun(dd, 2)));
  foreach([2,3,5,6,7,10,15], mm,
          addc(Str("log(", mm, ")^2"), log(mm)^2));
  foreach([2,3,5,6,7,10,15], mm,
          addc(Str("Pi*log(", mm, ")"), Pi*log(mm)));
  foreach([2,3,5,6,8,12], mm,
          addc(Str("Pi*log((sqrt(", mm, ")+1)/(sqrt(", mm, ")-1))"),
               Pi*log((sqrt(mm)+1)/abs(sqrt(mm)-1))));
  foreach([2,3,5,6,8,12], mm,
          addc(Str("log((sqrt(", mm, ")+1)/(sqrt(", mm, ")-1))^2"),
               log((sqrt(mm)+1)/abs(sqrt(mm)-1))^2));
  foreach([1,2,3,4,5,6,7,8], mm, addc(Str("Pi*atan(1/", mm, ")"), Pi*atan(1/mm)));
  foreach([1,2,3,4,5,6,7,8], mm, addc(Str("atan(1/", mm, ")^2"), atan(1/mm)^2));
  foreach([2,3,5,7], m1, foreach([2,3,5,7], m2,
      if(m1 <= m2, addc(Str("log(", m1, ")*log(", m2, ")"), log(m1)*log(m2)))));
  foreach([1/2,1/3,2/3,1/4,3/4,1/5,2/5,3/5,4/5,1/6,5/6,1/8,1/9,-1,-1/2,-1/3,-2,-3,1/16,1/27], zz,
      addc(Str("Li2(", zz, ")"), -intnum(tt = 0, 1, log(1-zz*tt)/tt)));
  addc("Gamma(1/4)^4/Pi", gamma(1/4)^4/Pi);
  addc("Gamma(1/3)^6/Pi^2", gamma(1/3)^6/Pi^2);
  foreach([2,3,5,6,7], mm, addc(Str("zeta(2)/sqrt(", mm, ")"), zeta(2)/sqrt(mm)));
  foreach([2,3,5], mm, addc(Str("Pi*log(2)/sqrt(", mm, ")"), Pi*log(2)/sqrt(mm)));
  foreach([2,3,5], mm, addc(Str("Pi*log(3)/sqrt(", mm, ")"), Pi*log(3)/sqrt(mm)));
}

HMAXR = 10^9;

try1(xv, i) =
{ my(rel, hh);
  rel = lindep([1, xv, VALS[i]]);
  if(rel == 0, return(0));
  if(rel[3] == 0 || rel[2] == 0, return(0));
  hh = vecmax(apply(abs, rel));
  if(hh > HMAXR, return(0));
  if(abs(rel[1] + rel[2]*xv + rel[3]*VALS[i]) > 10.0^(-PREC+10)*hh, return(0));
  [NAMES[i], rel];
}

try2(xv, i, j) =
{ my(rel, hh);
  rel = lindep([1, xv, VALS[i], VALS[j]]);
  if(rel == 0, return(0));
  if(rel[2] == 0, return(0));
  if(rel[3] == 0 && rel[4] == 0, return(0));
  hh = vecmax(apply(abs, rel));
  if(hh > 10^5, return(0));
  if(abs(rel[1] + rel[2]*xv + rel[3]*VALS[i] + rel[4]*VALS[j]) > 10.0^(-PREC+12)*hh, return(0));
  [Str(NAMES[i], " , ", NAMES[j]), rel];
}

identify(lab, xv) =
{ my(found = 0, r);
  r = lindep([1, xv]);
  if(r != 0 && r[2] != 0 && vecmax(apply(abs, r)) < 10^12
     && abs(r[1] + r[2]*xv) < 10.0^(-PREC+8),
     print(lab, "  RATIONAL  xi = ", -r[1], "/", r[2]); return(0));
  for(i = 1, #VALS,
    r = try1(xv, i);
    if(r != 0, print(lab, "  ", r[1], "   rel = ", r[2]); found = 1));
  if(found == 0,
    for(i = 1, #VALS, for(j = i+1, #VALS,
      r = try2(xv, i, j);
      if(r != 0, print(lab, "  ", r[1], "   rel = ", r[2]); found = 1))));
  if(found == 0, print(lab, "  UNIDENTIFIED"));
  found;
}

{
  my(ln = readstr("out/xi.txt"));
  for(i = 1, #ln,
    my(v = strsplit(ln[i], " "));
    if(#v >= 2, identify(v[1], eval(Str(v[2], "*1.0")))));
}
quit;
