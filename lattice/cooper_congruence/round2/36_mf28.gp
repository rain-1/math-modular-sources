\\ 36_mf28.gp -- Task (4): what can PARI's mf package do at level 28, weight 5/2?
default(parisize, 8000000000);
default(timer,1);
print("=== mfinit([28,5/2]) ===");
{ E = 0;
  iferr(E = mfinit([28,5/2]), err, print("  ERROR: ", err); E = 0);
}
if(type(E)=="t_VEC", print("  mfdim = ", mfdim(E)); print("  mfparams = ", mfparams(E)), print("  no space"));
print("");
print("=== mfinit([28,5/2],1) (cuspidal) and mfinit([28,5/2],4) (new) ===");
{ iferr(print("  Scusp dim = ", mfdim(mfinit([28,5/2],1))), err, print("  ERROR cusp: ", err)); }
{ iferr(print("  Sold dim  = ", mfdim(mfinit([28,5/2],3))), err, print("  ERROR: ", err)); }
print("");
print("=== dimensions of the holomorphic spaces M_{5/2+12r}(Gamma_0(28)) ===");
{ for(r=0,2, my(w=5/2+12*r);
   iferr(print("  r=",r," weight ",w,"  dim M = ", mfdim(mfinit([28,w])), "   dim S = ", mfdim(mfinit([28,w],1))),
         err, print("  r=",r," weight ",w," ERROR: ", err)));
}
print("");
print("=== does PARI support half-integral weight at all? sanity at level 4 ===");
{ iferr(my(M=mfinit([4,5/2])); print("  level 4 wt 5/2 dim M = ", mfdim(M), "  basis q-exp:");
        my(B=mfbasis(M)); for(i=1,#B, print("    ", mfcoefs(B[i],12))),
     err, print("  ERROR: ", err)); }
{ iferr(my(M=mfinit([4,1/2])); print("  level 4 wt 1/2 dim = ", mfdim(M)), err, print("  wt 1/2 ERROR: ", err)); }
print("");
print("=== Kohnen plus space support? ===");
{ iferr(my(M=mfinit([4,5/2],[1,1])); print("  plus-space flag accepted: dim=", mfdim(M)), err, print("  no plus-space flag: ", err)); }
{ print("  mfkohnenbasis available? "); iferr(my(M=mfinit([4,5/2],1)); print("   ", mfkohnenbasis(M)), err, print("   ERROR: ", err)); }
{ iferr(my(M=mfinit([28,5/2],1)); print("  mfkohnenbasis(level 28, cusp): ", mfkohnenbasis(M)), err, print("  ERROR: ", err)); }
quit;
