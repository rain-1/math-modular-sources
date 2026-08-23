mz0 = [-1,2,-1,-1; 1,-4,2,2; 1,8,-4,-3; -6,-19,9,7];
mzi = [-8,-5,2,1; 14,9,-4,-2; -22,-15,5,3; 45,32,-12,-7];
mz1 = [2,1,0,0; -3,-2,0,0; 0,0,1,0; 0,0,-3,-1];
jz  = [0,0,0,1; 0,0,1,1; 0,-1,0,0; -1,-1,0,0];
id4 = matid(4);
n0 = mz0^3-id4; ninf = mzi^6-id4;
s0 = jz*n0;

print("=== C4: induced rank-2 form on Z^4/rad(S0) ===");
kk = matkerint(s0);
sv = matsnf(kk,1); uu = sv[1]; dd = sv[3];
print("U*K*V = D = ", dd, "  (kernel maps to e3,e4)   det U = ", matdet(uu));
ui = uu^(-1);
gp2 = mattranspose(ui)*s0*ui;
print("full Gram in new basis = ", gp2);
gg = [gp2[1,1],gp2[1,2]; gp2[2,1],gp2[2,2]];
print("induced 2x2 form G = ", gg, "  det = ", matdet(gg));
a=gg[1,1]; b=2*gg[1,2]; c=gg[2,2];
print("as binary form (a,b,c) = (",a,",",b,",",c,")  disc = ", b^2-4*a*c, "  content = ", gcd([a,b,c]));
print("positive definite? a>0 and det>0 : ", a>0 && matdet(gg)>0);
print("reduced(G/content) = ", qfbred(Qfb(a/2,b/2,c/2)));
print("reduced(A2/content)= ", qfbred(Qfb(1,-1,1)), "   [A2=(2,-2,2)]");
print("qfminim(G,4) = ", qfminim(gg,4));
print("qfminim(A2,4)= ", qfminim([2,-1;-1,2],4));
tt = qflllgram(gg);
print("LLL Gram of G = ", mattranspose(tt)*gg*tt, "  (transform det ", matdet(tt), ")");
print("G ~ A2 over GL2(Z)? ", mattranspose(tt)*gg*tt == [2,-1;-1,2] || mattranspose(tt)*gg*tt == [2,1;1,2]);

print();
print("=== C7: 18 local kernel dimensions over Q(zeta_6) ===");
z6 = Mod('x, polcyclo(6));
print("check z6^6 = ", z6^6, "  z6^3 = ", z6^3, "  z6^2+... minpoly ok");
{
my(tot=0, tot2=0);
print("k | d0=dimker(mz0*z6^(2k)-I)  d1=dimker(mz1*z6^(3k)-I)  di=dimker(mzi*z6^(k)-I) | H^1 | di' (mzi^-1) -> H^1'");
for(k=0,5,
  my(a0 = 4 - matrank(mz0*z6^(2*k) - id4));
  my(a1 = 4 - matrank(mz1*z6^(3*k) - id4));
  my(ai = 4 - matrank(mzi*z6^(k)   - id4));
  my(ai2= 4 - matrank(mzi^(-1)*z6^(k) - id4));
  my(h = 4 - (a0+a1+ai)); my(h2 = 4-(a0+a1+ai2));
  tot += h; tot2 += h2;
  print(k, " |  ", a0, "  ", a1, "  ", ai, "  | ", h, " | ", ai2, " -> ", h2);
);
print("total dim H_par = ", tot, "   (alt convention: ", tot2, ")");
}
quit();
