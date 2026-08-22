\\ calegari_rows.gp -- reconstruct the rows of Calegari, "Irrationality of certain
\\ p-adic periods for small p" (math/0408214), sec.3 (zeta_p(3), X_0(p) genus 0)
\\ and sec.4 (L_2(2,chi_{-4}) = zeta_2(2), X_1(4)), and MEASURE
\\    sigma_p = lim v_p(a_n - b_n eta)/n            (p-adic decay of the linear form)
\\    log lambda_1 = lim log max(|a_n|,|b_n|)/n     (archimedean growth)
\\    k = smallest exponent with lcm(1..n)^k a_n in Z
\\ and check them against Calegari's asserted 12/(p-1), 6/(p-1), 2k+1.
\\ Run: gp -q /home/ubuntu/code/math-modular-sources/lattice/padic_irrationality/calegari_rows.gp
default(parisizemax,"6G");
default(realprecision,60);

NQ = 200;

etaprod(m,N) = { my(s=1+O(q^N)); for(n=1,(N-1)\m, s*=(1-q^(m*n))); s; }

\\ f = (Delta(p tau)/Delta(tau))^{1/(p-1)},  local parameter at i*oo on X_0(p)
unif(pp,N) = { q*(etaprod(pp,N)/etaprod(1,N))^(24/(pp-1)) + O(q^(N+1)); }

\\ E_2^*(tau) = E_2(tau) - p E_2(p tau),  E_2 = -1/24 + sum sigma_1(n) q^n
E2star(pp,N) = { my(s = (pp-1)/24 + O(q^(N+1)));
  for(n=1,N, s += (sigma(n,1) - if(n%pp==0, pp*sigma(n/pp,1), 0))*q^n); s; }

\\ E'_{-2} = theta^{-3}(E_4 - V_p E_4)
Eprim(pp,N) = { my(s = O(q^(N+1)));
  for(n=1,N, s += ((sigma(n,3) - if(n%pp==0, sigma(n/pp,3), 0))/n^3)*q^n); s; }

\\ H = E_2^*(E'_{-2} + eta) = sum (a_n + eta b_n) f^n   (Calegari's normalisation, x24)
zetarow(pp,N) = {
  my(f=unif(pp,N), qf, Es=E2star(pp,N), Ep=Eprim(pp,N), HA, HB, a, b);
  qf = serreverse(f);
  HA = subst(Es*Ep, q, qf); HB = subst(Es, q, qf);
  a = vector(N+1, j, 24*polcoef(HA,j-1));
  b = vector(N+1, j, 24*polcoef(HB,j-1));
  [a,b];
}

\\ ---- Calegari sec.4: X_1(4), z = (Delta(4 tau)/Delta(tau))^{1/3}, F_1, F'_{-1}
catrow(N) = {
  my(z, qz, F1, Fm, HA, HB, a, b);
  z = q*(etaprod(4,N)/etaprod(1,N))^8 + O(q^(N+1));
  \\ F_1 = L(0,chi)/2 + sum_{n>=1} (sum_{d|n} chi_{-4}(d)) q^n ; L(0,chi_{-4}) = 1/2
  F1 = 1/4 + O(q^(N+1));
  for(n=1,N, my(c=0); fordiv(n,d, c += kronecker(-4,d)); F1 += c*q^n);
  \\ F'_{-1} = sum_{n>=0} q^{2n+1}(-1)^n /((2n+1)^2 (1-q^{2n+1}))
  Fm = O(q^(N+1));
  for(m=0,(N-1)\2, my(r=2*m+1, s=O(q^(N+1))); for(j=1,N\r, s += q^(r*j));
      Fm += (-1)^m/r^2*s);
  qz = serreverse(z);
  HA = subst(F1*Fm, q, qz); HB = subst(F1, q, qz);
  a = vector(N+1, j, 4*polcoef(HA,j-1));
  b = vector(N+1, j, 4*polcoef(HB,j-1));
  [a,b];
}

padicdata(a,b,pp,N,nlist) = {
  my(x = a[N+1]/b[N+1], out=vector(#nlist));
  for(j=1,#nlist, my(n=nlist[j], v=valuation(a[n+1]-b[n+1]*x, pp));
      out[j] = [n, v, v*1./n]);
  out;
}
archdata(a,b,nlist) = {
  vector(#nlist, j, my(n=nlist[j]); [n, log(max(abs(a[n+1]*1.),abs(b[n+1]*1.)))/n]);
}
kexp(a,N) = {
  my(kmax=0,d=1);
  for(n=1,N, d=lcm(d,n); my(kk=0); while(denominator(a[n+1]*d^kk)!=1 && kk<12, kk++);
      if(kk>kmax,kmax=kk));
  kmax;
}

print("=== Calegari sec.3: zeta_p(3) rows on X_0(p), p = 2,3,5,7,13 ===");
{
for(i=1,5,
  my(pp=[2,3,5,7,13][i], N=NQ, ab, a, b, nl);
  ab = zetarow(pp,N); a=ab[1]; b=ab[2];
  nl = [N\4, N\2, 3*N\4];
  print("\np = ", pp);
  print("  b_0..b_6 = ", vector(7,j,b[j]));
  print("  a_0..a_6 = ", vector(7,j,a[j]));
  print("  k (lcm exponent for a_n, n<=60) = ", kexp(a,60), "   [Calegari: 2k+1 = 3]");
  print("  measured v_p(a_n - b_n eta)/n  : ", padicdata(a,b,pp,N,nl));
  print("     [Calegari sigma_p = 12/(p-1) = ", 12/(pp-1), "]");
  print("  measured log max(|a|,|b|)/n    : ", archdata(a,b,nl));
  print("     [Calegari log lambda_1 = (6/(p-1))log p = ", 6/(pp-1)*log(pp*1.), "]");
);
}

print("\n=== Calegari sec.4: the 2-adic Catalan row on X_1(4) ===");
{
  my(N=NQ, ab=catrow(N), a=ab[1], b=ab[2], nl=[N\4,N\2,3*N\4]);
  print("  b_0..b_6 = ", vector(7,j,b[j]), "   [paper: -1,-4,28,-272,3036,-36624,464368]");
  print("  a_0..a_6 = ", vector(7,j,a[j]), "   [paper: 0,1,-3,116/9,-331/9,-99116/225,3133076/225]");
  print("  k (lcm exponent, n<=60) = ", kexp(a,60), "   [Calegari: 2]");
  print("  measured v_2(a_n - b_n eta)/n : ", padicdata(a,b,2,N,nl), "   [Calegari: 8]");
  print("  measured log max(|a|,|b|)/n   : ", archdata(a,b,nl), "   [Calegari: 4 log 2 = ", 4*log(2.),"]");
}
quit;
