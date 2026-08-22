/* lattice/positivity/pairs.gp  --- Task 3 (b): TWO-ROW positive constructions
   over the moment grid.  For a pair S = {(m1,j1),(m2,j2)} of Catalan-world
   Beukers kernels put d = lcm of the four denominators, a_i = d*A_i,
   b_i = d*B_i in Z, h = a1*b2 - a2*b1 the mixed minor, and
        MOD = d * 2^{v_2(h)}     (the 2-adic bridge of CATALAN_AUDIT sec.4(d))
   or   MOD = d                  (canonical, no bridge).
   L_S = {p in Z^2 : MOD | p.a and MOD | p.b}; q = p.a/MOD, p' = p.b/MOD in Z.
   Both moments are POSITIVE, so on L_S cap Z^2_{>=0} the form
   q*G - p' = p1*M1 + p2*M2 > 0 is non-vanishing BY A THEOREM.
   We report the cone minimum, the unrestricted first minimum, and the ratio.
   Prepend rows_pos.gp and cone80.gp (for conemin/redu).                     */
\p 700

/* LLL reduction with a modest, adaptive scale (the 600-digit offset used in
   cone80.gp overflows here because MOD is comparatively small). */
{
redu2(BB, lz, ln) = my(sc = 10^(min(150, min(precision(1.*lz),precision(1.*ln))-40-ceil(log(1.*vecmax([abs(BB[1,1]),abs(BB[1,2]),abs(BB[2,1]),abs(BB[2,2])]))/log(10))) + ceil(log(1/min(lz,ln))/log(10))));
 BB*qflll([round(BB[1,1]*lz*sc), round(BB[1,2]*lz*sc);
           round(BB[2,1]*ln*sc), round(BB[2,2]*ln*sc)]);
}
{
pairstat(A1,B1,A2,B2, M1, M2, bridge, tag) =
 localprec(700);
 my(dd=lcm(lcm(denominator(A1),denominator(B1)),
           lcm(denominator(A2),denominator(B2))),
    a1=dd*A1, b1=dd*B1, a2=dd*A2, b2=dd*B2,
    hh=a1*b2-a2*b1, v2h=if(hh==0, 0, valuation(hh,2)),
    MOD=if(bridge, dd*2^v2h, dd));
 my(B0=klat2(a1,b1,a2,b2,MOD), idx=abs(matdet(B0)));
 /* for an integer coefficient vector p the form is
    (p.a*G - p.b)/MOD = d*(p1*M1 + p2*M2)/MOD, so the scaled weights are: */
 my(lz=dd*M1/MOD, ln=dd*M2/MOD);   /* both > 0 : cone = first quadrant */
 my(Br=redu2(B0, lz, ln));
 my(nrm(v)=sqrt((v[1]*lz)^2+(v[2]*ln)^2));
 my(v1=[Br[1,1],Br[2,1]], v2=[Br[1,2],Br[2,2]]);
 my(l1=nrm(v1), l2=nrm(v2));
 if(l1>l2, my(tm=l1); l1=l2; l2=tm);
 my(R=200, r1=conemin(Br,lz,ln,R), r2=conemin(Br,lz,ln,2*R), guard=0);
 while(guard<4 && (r1[1]==0 || r1[1]!=r2[1]), R*=2; r1=r2; r2=conemin(Br,lz,ln,2*R); guard++);
 my(fc=r2[1], cv=r2[2]);
 if(fc==0, return([0,0,0,0,0]));
 my(q=(cv[1]*a1+cv[2]*a2)/MOD, pp=(cv[1]*b1+cv[2]*b2)/MOD);
 my(ok = (type(q)=="t_INT") && (type(pp)=="t_INT")
        && abs(abs(q*Catalan-pp)-fc) < fc*10.^(-precision(1.*lz)+50));
 [fc, l1, fc/l1, if(q==0, 0, log(abs(q))), if(ok && q!=0,1,0), v2h, log(idx), log(MOD), log(dd)];
}
