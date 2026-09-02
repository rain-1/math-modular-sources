/* 02_poles.gp -- TASK 2.  Where Phi is meromorphic; exact pole locations, order of
   ramification nu, exact principal part, and the sharp asymptotics of c(m).      */
default(parisize, 1000000000);
read("lib.gp");
\p 60
NT = 500;
u7(t)  = (eta(7*t,1)/eta(t,1))^4;
u10(t) = (eta(5*t,1)*eta(10*t,1)/(eta(t,1)*eta(2*t,1)))^2;
u18(t) = (eta(2*t,1)*eta(3*t,1)^2*eta(18*t,1)/(eta(t,1)*eta(6*t,1)^2*eta(9*t,1)))^6;
UU = [u7,u10,u18];
DS = [[1,7],[1,2,5,10],[1,2,3,6,9,18]];
RS = [[-4,4],[-2,-2,2,2],[-6,6,12,-12,-6,6]];
T0 = [ (5+I*sqrt(3))/14, (3+I)/10, (3+I)/6 ];
U0 = [ (-13+3*sqrt(-3))/98, (-3+4*I)/25, -7+4*sqrt(3) ];
NAM = ["s7","s10","s18"];
{ E2n(t) = my(q=exp(2*Pi*I*t), s=1, p=1); for(n=1,NT, p*=q; s -= 24*sigma(n)*p); s; }
{ Fn(k,t) = my(s=0,d=DS[k],r=RS[k]); for(j=1,#d, s += r[j]*d[j]*E2n(d[j]*t)); s/24; }
{ Phin(k,t) = my(u=UU[k](t), F=Fn(k,t), B=ROWS[k][3], C=ROWS[k][4]); F^2*u*(1-C*u^2)/(1+B*u+C*u^2)^2; }
{ ligozat(N, ds, rs) = my(v=[]); fordiv(N,c, my(s=0); for(j=1,#ds, s += gcd(c,ds[j])^2*rs[j]/(gcd(c,N\c)*c*ds[j])); v=concat(v,[[c,eulerphi(gcd(c,N\c)),N/24*s]])); v; }
print("=== divisor of u at the cusps  [c, #cusps, ord] ; degree of u on X_0(N) ===");
{ for(k=1,3, my(v=ligozat(ROWS[k][2],DS[k],RS[k]), s=0);
  for(j=1,#v, if(v[j][3]>0, s += v[j][2]*v[j][3]));
  print("  ",NAM[k],": ",v,"   deg u = ",s)); }
print();
print("=== roots of 1+Bu+Cu^2 and the points of H where u attains them ===");
{ for(k=1,3, my(g = 1+ROWS[k][3]*'u+ROWS[k][4]*'u^2);
  print("  ",NAM[k],"  roots = ",polroots(g));
  print("      u(tau_0) at tau_0 = ",T0[k]," : ", UU[k](T0[k]), "   1+Bu+Cu^2 = ", subst(g,'u,UU[k](T0[k])));
  print("      minimal poly of tau_0 : ", ROWS[k][2], "t^2 ...  disc = ", [-3,-4,-36][k])); }
print();
print("=== order nu = ord_{tau-tau_0}(u-u_0)  (numerical) ===");
{ for(k=1,3, my(r=[]);
  for(j=5,9, my(e=10.0^(-j)); r=concat(r,[round(1000*log(abs(UU[k](T0[k]+2*e)-U0[k])/abs(UU[k](T0[k]+e)-U0[k]))/log(2.0))/1000.]));
  print("  ",NAM[k],"  nu = ",r)); }
print();
print("=== principal part  Phi = A2 y^-2 + A1 y^-1 + O(1),  y = tau - tau_0 ===");
print("    predicted  A2 = nu^2/(4 pi^2 g'(u0)),  g'(u0)^2 = B^2-4C = lambda_1 lambda_2,  A1 = i A2 / Im(tau_0)");
NU = [3,4,4];
{ for(k=1,3, my(u0=U0[k], gp0=ROWS[k][3]+2*ROWS[k][4]*u0, pred=NU[k]^2/(4*Pi^2*gp0), e=10.0^(-8), a2, a1, y0=imag(T0[k]));
  a2 = e^2*Phin(k,T0[k]+e); a1 = e*(Phin(k,T0[k]+e)-pred/e^2);
  print("  ",NAM[k],"  g'(u0) = ",gp0,"   A2/pred = ",a2/pred,"   A1/A2 = ",a1/pred,"   i/Im(tau0) = ",I/y0)); }
print();
print("=== growth of c(m): R = 1/|q_0| = exp(2 pi Im tau_0) and the sharp asymptotic ===");
M = 400; NN = M+6;
{ CV = vector(3); for(k=1,3, my(S=Setup(k,NN)); CV[k]=vector(M,m,polcoeff(S[4],m))); }
{ for(k=1,3, my(t0=T0[k],u0=U0[k],gp0=ROWS[k][3]+2*ROWS[k][4]*u0,q0=exp(2*Pi*I*t0),y0=imag(t0),A=-NU[k]^2/gp0,r=[]);
  print("  ",NAM[k],"  R = ",1/abs(q0),"   |c(400)|^(1/400) = ",abs(CV[k][400]*1.0)^(1/400));
  for(j=0,3, my(m=396+j, pr=A*(m-1/(2*Pi*y0))*q0^(-m));
     pr = if(k<3, 2*real(pr), real(pr));
     if(abs(pr)>10^20, r=concat(r,[m, CV[k][m]/pr])));
  print("      c(m) / [ -nu^2 (q0^-m/g'(u0) + c.c.)(m - 1/(2 pi Im tau0)) ] : ", r)); }
quit;
