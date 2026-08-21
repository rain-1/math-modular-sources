/* 05_periodpoly.gp -- the Fricke period polynomial of the four interior-fold
   rows is  R(tau) = L(Phi,3) + 2 pi i L(Phi,2) tau - 2 pi^2 L(Phi,1) tau^2
   and the W_N-invariance of H_xi forces  R(tau) = xi (1 + N tau^2).
   Hence  L(Phi,2) = 0  and  L(Phi,1) = -N L(Phi,3)/(2 pi^2).              */
read("common.gp");
default(realprecision, 50);
NL = [["alpha",12],["gamma",6],["eps",8],["zeta",9]];
{
for(i=1,#NL,
  my(nm=NL[i][1], N=NL[i][2], L3=Ltarget(nm), L1);
  \\ L(Phi,1) = P(1)*lim_{s->1} L(psi,s)L(phi,s-3);  psi=phi=1 except zeta
  L1 = if(nm=="zeta", lfun(-3,1)*lfun(-3,-2), Pmellin(nm,1)*zetahurwitz(-2,1,1));
  print(nm,"  N=",N,"   L(Phi,1)=",L1,"   -N*L(Phi,3)/(2 Pi^2)=",-N*L3/(2*Pi^2),
        "   diff=", L1 + N*L3/(2*Pi^2), "    L(Phi,2)=", LPhi(nm,2));
);
}
quit;
