# Magnetic sources on Apéry's host: the source exists, keeps the free integration, and its pole sits on the unit circle

*Fable, 2026-09-02. Closure of top-9 item 1. Full report: `lattice/magnetic_apery/REPORT.md`.*

## 1. The idea and the census

A magnetic source (integral antiderivative, one free integration) on Apéry's $\zeta(3)$ host would give a depth-3 rate, $(\sqrt2-1)^4=0.0294$, with depth-2 denominators $d_n^2$: CDT's own budget on a far better host. On $X_0(6)$ the census is complete: $u$ is a Hauptmodul, $\Gamma_0(6)$ has no elliptic points, $\operatorname{div}F=(\tfrac12)+(\tfrac13)$, so $\Phi=F^2\rho(u)$ with $\rho\in\mathbf Q(u)$, $\rho(0)=0$, exhausts the weight-4 meromorphic forms with $a_0=0$. An exact lattice descent with an LLL detector over every CM orbit of discriminant $\ge-96$, all pairs and triples, pole orders $1$–$4$, finds **exactly three** families of magnetic forms.

## 2. The hit and the obstruction

**The source.** $\Phi_{-8}=F^2\cdot\dfrac{-u(1+18u+72u^2)(72u^2-1)}{(1+16u+72u^2)^2(1+17u+72u^2)}$, double pole at the discriminant-$(-8)$ CM point $\tau_0=-\tfrac13+i\sqrt2/6$, $W_6$-antiinvariant, $m\mid c(m)$ to $m=1000$. It is the canonical source of the *other* level-6 Fricke host, $(C,B)=(81,14)$. Transplanted to Apéry's host it keeps $k=2$ ($n\le250$) and $B_n/A_n\to\zeta(2)/8=\pi^2/48$ (185 digits).

**Theorem U (the unit obstruction).** $1/x=W+B$ with $W=Cu+1/u$ a Hauptmodul with expansion $1/q+\mathbf Z[[q]]$, so singular values of $W$ are algebraic integers and polar divisors are Galois-stable; hence $\max_i|W_i+B|\ge|N(W_0+B)|\ge1$, i.e. **the companion of a magnetic source has a singularity at $|x|\ge1$**, where $x(\tau_0)$ is the coordinate of the pole. On Apéry's host $N=\pm1$ on all 14 CM orbits and $x(\tau_0)=1$ exactly. A pole of $\Phi$ at $|x|=1$ puts a singularity of $B-\xi F$ at radius $1$, far inside $|t_2|=33.97$, so $A_n\xi-B_n\asymp n^{-2}$ and $d_n^2|A_n\xi-B_n|\to\infty$: no irrationality, whatever $k$ is. The only pole locus the coordinate cannot see is $x=\infty$, which on Apéry's host is the cusp pair, and there every $W_6=-1$ magnetic form is Bol-trivial ($\Phi=D^3g$, $B_n/A_n$ exactly constant, rational limits).

**Why Cooper's rows work.** For $s_7,s_{10},s_{18}$ the pole of the source is at $1+Bu+Cu^2=0$, i.e. at $x=\infty$ of the row's own coordinate: a magnetic source only helps the host it is native to. Verified on eight hosts: $N(W+B)=0$ exactly at the host's own pole discriminant ($-3$ for $s_7$, $-4$ for $s_{10}$, $-8$ for $C=81$, $-12$ for $C=64$, $-36$ for $s_{18}$, $-96$ for the level-12 host), nonzero integer elsewhere. Five further magnetic forms found on other hosts, all with $|W+B|\ge4$.

**Prize forgone.** A usable magnetic source would have scored entry $+0.846$, margin $+0.005$: CDT's own passing score on Apéry's host, against $-0.90$ / $-24$ for the real $k=3$ row; classically it would have given $\mu(\zeta(2))\le4.62$ against the record $5.09$.

## 3. By-products

* **Free integration $\iff$ meromorphy, 12 of 12 hosts.** $\limsup|c(m)|^{1/m}>1.87$ on the eight hosts with $k=2$ and $=1$ on the four without (Apéry, ε, ζ, Domb): Paşol–Zudilin's folklore conjecture and its converse, on this census.
* $L\bigl(D^{-1}(\Delta/E_4^2),2\bigr)=0.06103925100753799194274901793938221534\ldots$ is not a rational multiple of $1,\zeta(2),\zeta(3),\pi^2,\pi^3$ and not algebraic of degree $\le4$: an apparently new constant, the depth-2 period of Paşol–Zudilin's level-one magnetic form. For $E_4\Delta/E_6^2$ the pole sits at $\tau=i$ on the Fricke geodesic and the period integral diverges.

## 4. Verdict

Arithmetically right, geometrically dead: magnetic sources have their poles at CM points, the Apéry coordinate at a CM point is an algebraic unit, and a pole on the unit circle swamps any exponential rate. Item 1 is closed.
