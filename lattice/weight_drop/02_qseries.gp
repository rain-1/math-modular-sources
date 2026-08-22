/* Cooper s7,s10,s18: q-series of x,t,F,Phi,Xi=theta^{-1}Phi, g=sqrt(F), Psi_root=g^3 u.
   Verifies: F=sum A_n t^n ; Phi=F*theta_q t ; Xi integral (Bogner (n+1)|A_n) ;
             Theta_par = theta^{-2} Xi  and  B_n=[t^n](F*Theta_par).            */
NQ = 60;
ee(d) = eta(q^d + O(q^NQ));
E2s(d) = 1 - 24*sum(n=1,(NQ-1)\d, sigma(n)*q^(d*n)) + O(q^NQ);
th(f) = q*deriv(f,q);
{invth(f,k) = my(M=NQ-2); sum(m=1,M, polcoeff(f,m)/m^k*q^m) + O(q^(M+1));}
/* t-expansion of a q-series Y in powers of T (T = q + ...) */
{tco(Y,T,nmax) = my(co=vector(nmax+1),R=Y); for(n=0,nmax, co[n+1]=polcoeff(R,0); R=R-co[n+1]; if(n<nmax,R=R/T)); co;}

XX7 = (ee(1)/ee(7))^4;            T7  = q*XX7/(XX7^2+13*q*XX7+49*q^2);
F7  = (1/6)*(7*E2s(7)-E2s(1));
XX10= (ee(2)*ee(5)/(ee(1)*ee(10)))^6; T10 = q*XX10/(XX10-q)^2;
F10 = (1/12)*(10*E2s(10)+5*E2s(5)-2*E2s(2)-E2s(1));
XX18= ee(3)^4*ee(6)^4/(ee(1)^2*ee(2)^2*ee(9)^2*ee(18)^2); T18 = q*XX18/(XX18+3*q)^2;
F18 = (1/4)*(18*E2s(18)-9*E2s(9)-12*E2s(6)+6*E2s(3)+2*E2s(2)-E2s(1));

{par(ab,N) = my(a=ab[1],b=ab[2],c=ab[3],d=ab[4],A=vector(N+1),B=vector(N+1));
  A[1]=1;A[2]=b;B[1]=0;B[2]=1;
  for(n=1,N-1, A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]+n*(c*n^2+d)*A[n])/(n+1)^3;
               B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]+n*(c*n^2+d)*B[n])/(n+1)^3;);[A,B];}

dat = [["s7",T7,F7,[13,4,27,-3],7,1],["s10",T10,F10,[6,2,64,-4],10,2],["s18",T18,F18,[14,6,-192,12],18,2]];
NT = 22;
{for(r=1,#dat,
  my(nm=dat[r][1],T=dat[r][2],F=dat[r][3],N=dat[r][5],lam=dat[r][6]);
  my(PB=par(dat[r][4],NT+2), A=PB[1], B=PB[2]);
  print("======== ",nm,"   N=",N);
  print(" t = ",T+O(q^8));
  print(" F = ",F+O(q^8));
  my(cF=tco(F,T,NT));
  print(" F = sum A_n t^n ?  max diff = ", vecmax(apply(abs,vector(NT+1,k,cF[k]-A[k]))));
  my(PHI=F*th(T));  print(" Phi = ",PHI+O(q^9));
  my(XI=invth(PHI,1));
  print(" Xi = theta^{-1}Phi = ",XI+O(q^10));
  print(" Xi integral?  denominators: ", vector(12,m,denominator(polcoeff(XI,m))));
  my(THE=invth(PHI,3), TH2=invth(XI,2));
  print(" theta^{-3}Phi = theta^{-2}Xi ?  ", THE-TH2);
  my(cB=tco(F*THE,T,NT));
  print(" B_n = [t^n](F theta^{-2}Xi) ?  max diff = ", vecmax(apply(abs,vector(NT+1,k,cB[k]-B[k]))));
  my(g=sqrt(F));
  print(" g=sqrt(F) = ",g+O(q^12));
  print(" g denominators: ", vector(14,m,denominator(polcoeff(g,m))));
  my(PSI=g^3*T/lam);
  print(" Psi_root=g^3 t/lam = ",PSI+O(q^12));
  print(" Psi denominators: ", vector(14,m,denominator(polcoeff(PSI,m))));
);}
quit;
