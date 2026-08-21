/* G_2 (=8*zeta_2(2) normalisation of the repo) from the Zudilin row, 2-adically. */
{
G2rat(M) = my(r=zud(M)); [r[2][M+1]/r[1][M+1], 8*M-1-4*hammingweight(M)];
}
/* v_2 of a rational */
v2(x) = if(x==0, +oo, valuation(x,2));
