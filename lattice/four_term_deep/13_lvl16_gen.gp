/* Regenerate the level-16 Catalan host A(x) and its Phi_1 = E companion B(x). */
default(parisize, 3000000000);
if(type(NTERM)!="t_INT", NTERM = 200);
read("/home/ubuntu/code/math-modular-sources/lattice/catalan_two_classes/00_setup.gp");
H16 = mkhost(x16,F16);
Av  = Aof(H16);
Bv  = Bof(H16, mkPhi(Sin,[[1,1]]));   /* Phi_1 = E */
NMAX = NT-2;
print("NTERM = ", NTERM, "   NMAX = ", NMAX);
print("A first 15: ", vector(15,i,polcoeff(Av,i-1)));
Aref = [1,0,-4,16,-44,96,-176,320,-860,3520];
print("A[0..9] matches quoted A(x)? ", vector(10,i,polcoeff(Av,i-1))==Aref);
print("B first 15: ", vector(15,i,polcoeff(Bv,i-1)));
/* integrality of a_n */
{my(bad=0); for(n=0,NMAX, if(denominator(polcoeff(Av,n))!=1, bad=n; break)); 
 print("a_n integral for n<=",NMAX,"? ", bad==0, if(bad,Str(" first non-integer n=",bad),""));}
write1("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_A.txt","");
{for(n=0,NMAX, write("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_A.txt", polcoeff(Av,n)));}
write1("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_B.txt","");
{for(n=0,NMAX, write("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/out/lvl16_B.txt", polcoeff(Bv,n)));}
print("wrote out/lvl16_A.txt and out/lvl16_B.txt");
default(realprecision,60);
print("b_n/a_n at n=NMAX: ", polcoeff(Bv,NMAX)/polcoeff(Av,NMAX)*1.0);
print("-G/2            : ", -Catalan/2);
