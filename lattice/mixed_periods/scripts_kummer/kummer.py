"""Kummer-type hosts H = (1 -+ N x)^{-1/k}: fold periods, closed forms.

Family 'M' (minus):  H = (1 - N x)^{-1/k}, fold delta = +1/N, D = N-1,
   Q(u) = u^k + D  (sigma=+1),  alpha_j = D^{1/k} e^{i pi (2j+1)/k}
Family 'P' (plus):   H = (1 + N x)^{-1/k}, fold delta = -1/N, D = N+1,
   Q(u) = u^k - D  (sigma=-1),  alpha_j = D^{1/k} e^{2 pi i j/k}
In both cases  1-t = |Q(u)|/N  and  c^{(a)}[kappa] = k int_0^1 u^{a-1} kappa /Q(u) du
for kappa = 1/(1-t) [i.e. N/|Q|] ... see period() below.
"""
from mpmath import mp, mpf, mpc, log, exp, pi, quad, polylog, im, re, sqrt, atan, fabs, mpmathify

def setup(dps):
    mp.dps = dps

class Host:
    def __init__(self, k, N, fam):
        self.k = k; self.N = N; self.fam = fam
        self.sigma = 1 if fam == 'M' else -1
        self.D = N - 1 if fam == 'M' else N + 1
        self.delta = mpf(1)/N if fam == 'M' else -mpf(1)/N
    def alphas(self):
        k, D = self.k, self.D
        r = mpf(D)**(mpf(1)/k)
        if self.fam == 'M':
            return [r*exp(mpc(0,1)*pi*(2*j+1)/k) for j in range(k)]
        else:
            return [r*exp(2*mpc(0,1)*pi*j/k) for j in range(k)]
    def Q(self, u):
        return u**self.k + self.sigma*self.D
    def absQ(self, u):
        return abs(self.Q(u))      # = N*(1-t)
    # ---- fold periods by u-quadrature ----
    def per_u(self, a, kern):
        k, N, D = self.k, self.N, self.D
        s = self.sigma
        if kern == 'B':      # kappa = 1/(1-t)
            f = lambda u: u**(a-1) / (u**k + s*D) * k
        elif kern == 'D':    # kappa = log(1-t)/(1-t)
            f = lambda u: u**(a-1) * log(abs(u**k+s*D)/N) /(u**k+s*D) * k
        elif kern == 'L':    # kappa = log(1-t)
            f = lambda u: u**(a-1) * log(abs(u**k+s*D)/N)
        elif kern == 'tB':   # kappa = t/(1-t)
            f = lambda u: u**(a-1) * (1 - abs(u**k+s*D)/N) /(u**k+s*D) * k
        elif kern == 'tD':   # kappa = t log(1-t)/(1-t)
            f = lambda u: u**(a-1) * (1-abs(u**k+s*D)/N) * log(abs(u**k+s*D)/N) /(u**k+s*D)*k
        else: raise ValueError(kern)
        v = quad(f, [0,1], maxdegree=10)
        if kern in ('L',):
            v = v * self.sigma / (self.N/mpf(k))   # see derivation: dt = sgn u^{k-1}du/(N/k)
        return v
    # ---- fold periods by direct t-quadrature (independent check) ----
    def per_t(self, a, kern):
        k, N = self.k, self.N
        sg = 1 if self.fam=='M' else -1
        H = lambda t: (1 - sg*N*t)**(-mpf(k-a)/k)
        if kern == 'B':   g = lambda t: H(t)/(1-t)
        elif kern == 'D': g = lambda t: H(t)*log(1-t)/(1-t)
        elif kern == 'L': g = lambda t: H(t)*log(1-t)
        elif kern == 'tB':g = lambda t: H(t)*t/(1-t)
        elif kern == 'tD':g = lambda t: H(t)*t*log(1-t)/(1-t)
        else: raise ValueError(kern)
        return quad(g, [0, self.delta], maxdegree=12)
    # ---- closed forms ----
    def w(self):
        """winding: sum_l arg(u-alpha_l) = pi*w, constant on [0,1]"""
        al = self.alphas(); u = mpf(1)/2
        tot = sum(mp.arg(u-a) for a in al)
        return tot/pi
    def Lam(self, j):
        al = self.alphas(); a = al[j]
        return log(1-a) - log(-a)
    def cB(self, a):
        al = self.alphas(); D = self.D; s = self.sigma
        return (-s/mpf(D))*sum(al[j]**a * self.Lam(j) for j in range(self.k))
    def Ipair(self, j, l):
        """int_0^1 Log(u-alpha_l)/(u-alpha_j) du, principal branches"""
        al = self.alphas(); aj, alph = al[j], al[l]
        if j == l:
            return (log(1-alph)**2 - log(-alph)**2)/2
        d = aj - alph
        P1 = (1-alph)/d          # z(1)
        P0 = (-alph)/d           # z(0)
        return (log(1-alph)*log(1-P1) + polylog(2,P1)) - (log(-alph)*log(1-P0) + polylog(2,P0))
    def cD(self, a):
        al = self.alphas(); D = self.D; s = self.sigma; k = self.k
        S = sum(al[j]**a * self.Ipair(j,l) for j in range(k) for l in range(k))
        return (-s/mpf(D))*S - (log(self.N) + mpc(0,1)*pi*self.w())*self.cB(a)
    def cL(self, a):
        """closed form for kappa=log(1-t): elementary logs"""
        al = self.alphas(); k=self.k; N=self.N
        tot = mpc(0)
        for l in range(k):
            alph = al[l]
            # int_0^1 u^{a-1} Log(u-alpha) du
            inner = sum(alph**i/mpf(a-i) for i in range(a)) + alph**a*(log(1-alph)-log(-alph))
            tot += (log(1-alph) - inner)/mpf(a)
        tot += -(log(N) + mpc(0,1)*pi*self.w())/mpf(a)
        return self.sigma*mpf(k)/N * tot
