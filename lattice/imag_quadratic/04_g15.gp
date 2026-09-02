\\ 04_g15.gp -- Gamma_1(5) is the ONLY four-point genus-zero host carrying a nebentypus
\\ of order > 2 (03_chars.out): the odd quartic psi_4, values in Q(i).  Does the
\\ fold-regular subspace  V = ker a_0(oo) cap ker a_0(cusp 0)  meet the psi_4-component?
\\ If it did, Gamma_1(5) would carry a Q(i)-rational Apery system with the non-real
\\ period L(2,psi_4), scored exactly like CDT's theorem and with NO number-field tax.
default(realprecision,60);
cycord(m) = {my(n,v=variable(m)); for(n=1,300, if(polcyclo(n,v)==m, return(n))); 0;}
toC(z) = {if(type(z)!="t_POLMOD", return(z*1.0)); my(m=z.mod, n=cycord(m)); subst(lift(z), variable(m), exp(2*Pi*I/n));}
G = znstar(5,1);
print("characters mod 5 (Conrey label m): order, chi(-1) as an element of Q/Z");
for(m=1,4, print("  m=",m,"  order=",charorder(G,znconreylog(G,m)),"  chi(-1)=",chareval(G,znconreylog(G,m),-1)));

print();
print("--- the full weight-3 Eisenstein space of Gamma_1(5) (both nebentypus components) ---");
mf2 = mfinit([5,3,[G,znconreylog(G,2)]],4);
mf3 = mfinit([5,3,[G,znconreylog(G,3)]],4);
B2 = mfbasis(mf2);
B3 = mfbasis(mf3);
print("dim M_3^Eis(Gamma_1(5)) = dim(psi_4-comp) + dim(psi_4bar-comp) = ",#B2,"+",#B3," = ",#B2+#B3);
MA = matrix(2,#B2+#B3);
for(j=1,#B2, MA[1,j]=toC(mfcoefs(B2[j],0)[1]));
for(j=1,#B2, my(p=0); MA[2,j]=toC(mfslashexpansion(mf2,B2[j],[0,-1;1,0],0,1,&p)[1]));
for(j=1,#B3, MA[1,#B2+j]=toC(mfcoefs(B3[j],0)[1]));
for(j=1,#B3, my(p=0); MA[2,#B2+j]=toC(mfslashexpansion(mf3,B3[j],[0,-1;1,0],0,1,&p)[1]));
print("rank of (a_0 at oo, a_0 at cusp 0) on the full 4-dimensional space = ",matrank(MA));
print("dim(fold-regular) = ",#B2+#B3-matrank(MA),"    [hostscan REPORT sec.10.3: 2, spanned over Q(sqrt5) by Phi_D and Phi_new]");

print();
print("--- restricted to the psi_4-component (Conrey m=2, order 4, odd) ---");
CHI = znconreylog(G,2);
mf = mfinit([5,3,[G,CHI]],4);
B = mfbasis(mf);
print("dim M_3^Eis(Gamma_1(5), psi_4) = ",#B);
print("basis 1 q-expansion : ",mfcoefs(B[1],5));
print("basis 2 q-expansion : ",mfcoefs(B[2],5));
M = matrix(2,#B);
for(j=1,#B, M[1,j]=toC(mfcoefs(B[j],0)[1]));
for(j=1,#B, my(p=0); M[2,j]=toC(mfslashexpansion(mf,B[j],[0,-1;1,0],0,1,&p)[1]));
print("a_0 at oo      : ",vector(#B,j,M[1,j]));
print("a_0 at cusp 0  : ",vector(#B,j,M[2,j]));
print("rank = ",matrank(M),"   determinant = ",matdet(M));
print("dim(fold-regular cap psi_4-component) = ",#B-matrank(M));
print();
print("Reading: basis 1 is the INNER direction E_3^{psi_4,1} (a_0 = 0 at oo, nonzero at 0),");
print("basis 2 is the OUTER direction E_3^{1,psi_4} (a_0 nonzero at oo, zero at 0).");
print("The two conditions are therefore independent ON the two-dimensional component and");
print("kill it entirely: there is no Q(i)-rational fold-regular source on Gamma_1(5).");
print("The 2-dimensional fold-regular subspace of the FULL space is the 'diagonal' one,");
print("stable under coefficientwise conjugation, and its interesting line is Phi_new,");
print("defined over the REAL quadratic field Q(sqrt5).");
quit;
