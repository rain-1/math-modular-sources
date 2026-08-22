/* lattice/positivity/rows_pos.gp
   Shared exact-row layer for the positivity programme.
   Rebuilt from lattice/catalan_audit/{rows.gp,nest2.gp} and
   lattice/catalan_positivity/rows_common.gp, with the moment generaliser of
   lattice/catalan_explicit/moments.gp.  No builtin names are shadowed
   (avoid psi, M, Phi, S, cmp: local matrices are called MAT, moduli MOD). */

default(parisizemax, 12000000000);

ee(m) = min(6*m, 4*m+3+logint(2*m-1,2));

/* --- Zudilin's Catalan row Q_m, P_m  (Zudilin 2002, Thm 1) --------------- */
{
zud(TOP) = my(QQ=vector(TOP+1), PP=vector(TOP+1), aa,bb,cc);
  QQ[1]=1; QQ[2]=7/4; PP[1]=0; PP[2]=13/8;
  for(m=1,TOP-1,
    aa=(2*m+1)^2*(2*m+2)^2*(20*m^2-8*m+1);
    bb=3520*m^6+5632*m^5+2064*m^4-384*m^3-156*m^2+16*m+7;
    cc=(2*m-1)^2*(2*m)^2*(20*m^2+32*m+13);
    QQ[m+2]=(bb*QQ[m+1]+cc*QQ[m])/aa;
    PP[m+2]=(bb*PP[m+1]+cc*PP[m])/aa);
  [QQ,PP];
}

/* --- General Catalan-world Beukers moment ------------------------------- *
 * mom(m,j) = int_0^1 int_0^1 x^{m-1/2+j}(1-x)^m y^{m+j}(1-y)^{m-1/2}
 *                    (1-xy)^{-(m+1+j)} dx dy
 *          = I(a,b,a,b,a+1)   with a = m+j, b = m,   0 <= j <= m.
 * Returns [AA,BB] with mom = AA*G - BB,  AA,BB in Q.
 * (j=0: Zudilin index m.  (m,j)=(3n,n): Nesterenko (4,7) index n.)          */
{
mom(m,jj) =
 my(NN=m, TOP=m+jj,
    cc=m!*prod(i=0,m-2,3/2+i)/TOP!,
    num=cc*prod(k=0,TOP-1,('t-k)), tail=vector(NN+2), pre=vector(NN+2),
    MAT, rhs, sol);
 pre[1]=1; for(j=1,NN, pre[j+1]=pre[j]*('t+j));
 tail[NN+2]=1; forstep(j=NN,0,-1, tail[j+1]=tail[j+2]*('t+j+1/2)^2);
 my(bs=vector(2*(NN+1)), lin=vector(NN+2));
 lin[1]=1; for(j=0,NN, lin[j+2]=lin[j+1]*('t+j+1/2));
 for(j=0,NN, bs[2*j+1]=pre[j+1]*lin[j+2]*tail[j+2];
             bs[2*j+2]=pre[j+1]*lin[j+1]*tail[j+2]);
 my(dg=2*NN+1);
 MAT = matrix(dg+1, 2*(NN+1), r, c, polcoeff(bs[c], r-1, 't));
 rhs = vectorv(dg+1, r, polcoeff(num, r-1, 't));
 sol = matsolve(MAT, rhs);
 my(A1=vector(NN+1,j,sol[2*(j-1)+1]), A2=vector(NN+1,j,sol[2*(j-1)+2]));
 my(BB = sum(j=0,NN, A2[j+1]),
    CC = sum(j=1,NN, sum(v=0,j-1, (v!/prod(i=0,v-1,3/2+i))*(A1[j+1]+A2[j+1]/(v+1/2)))));
 [4*BB, CC];
}

/* --- The two integer source rows, Lean/NearCriticalAssembly normalisation - */
{
zudrow(n) = my(m=3*n, DD=lcm(vector(6*n,i,i)), ZZ=zud(m));
 [2^ee(m)*DD^2*ZZ[1][m+1], 2^ee(m)*DD^2*ZZ[2][m+1]];
}
{
nestrow(n) = my(DD=lcm(vector(6*n,i,i)), r=mom(3*n,n));
 /* mom(3n,n) = (V_n G - U_n)/(4^{7n} D_{6n}^2), see CATALAN_EXPLICIT sec.1  */
 [4^(7*n)*DD^2*r[1], 4^(7*n)*DD^2*r[2]];
}

/* --- congruence lattice K_n = {c : M | c.X-row and M | c.Y-row} ---------- */
{
klat2(XX,YY,VV,UU,MOD) = my(K=matkerint([XX,VV,MOD,0,0,0; YY,UU,0,MOD,0,0]));
  mathnf(matrix(2,#K[1,],i,j,K[i,j]));
}
