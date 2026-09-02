P=Pol([2985984,-99035136,-80870400,165680640,306913536,216563328,82569024,18046944,2131344,95880,-3900,-398,1],'v);
Q=Pol([2985984,11446272,19761408,20307456,13828608,6571776,2234464,547648,96032,11752,953,46,1],'v);
print("P = ",P);
print("Q = ",Q);
print("factor(P) = ",factor(P));
print("factor(Q) = ",factor(Q));
print("gcd(P,Q) = ",gcd(P,Q));
print("content P = ",content(P),"  content Q = ",content(Q));
\\ evaluate at v_c = 1/(2 sqrt3)
default(realprecision,60);
vc=1/(2*sqrt(3));
print("R(v_c) = ",subst(P,'v,vc)/subst(Q,'v,vc));
print("U(tau_c) = R(v_c)*(Dv)^3/v_c^3 with Dv=0.4828903625506029613196821128714583372784332839665252579188:");
Dv=0.4828903625506029613196821128714583372784332839665252579188;
print("  = ",subst(P,'v,vc)/subst(Q,'v,vc)*Dv^3/vc^3);
print("  (should be U(tau_c)=128.8501310795959555981613385209584807459051702400753422862)");
w=Mod('t,'t^2-3);
print("R(v_c) exact in Q(sqrt3)? algdep: ",algdep(subst(P,'v,vc)/subst(Q,'v,vc),2));
quit;
