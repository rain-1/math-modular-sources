# A Non-Eisenstein Catalan Programme from the Eskandari–Murty–Nemoto Motive

## Abstract

Catalan's constant
\[
G=\beta(2)=\sum_{n\ge0}\frac{(-1)^n}{(2n+1)^2}
\]
is not known to be irrational. Recent work of Eskandari, Murty and Nemoto constructs a two-dimensional mixed motive over \(\mathbf Q\) having \(G\) as a period and produces systematic linear forms in \(1\) and \(G\).

This note records a series of calculations motivated by combining that construction with Apéry-style recurrence methods, Padé reduction, arithmetic holonomy and \(p\)-adic realization techniques. Several new structural phenomena emerge.

The principal results are:

1. an explicit low-order Picard–Fuchs equation for a natural moving EMN period;
2. an exact cubic distribution law explaining the factor \(10/9\) attached to the conductor-\(12\) Catalan realization;
3. a new non-Eisenstein real/\(3\)-adic realization of the Catalan extension;
4. explicit higher-pole recurrences and Catalan coefficients;
5. several natural Apéry-type systems and Padé/Jacobi/Racah constructions;
6. a family of homogeneous saddle-killing integrals with a positive *primary* irrationality score;
7. a precise explanation of why these constructions nevertheless do not currently prove \(G\notin\mathbf Q\).

The outcome is not an irrationality proof. It is, however, a rather concrete new research programme in which the remaining obstruction can be identified as arithmetic rather than geometric.

---

## 1. The geometric starting point

Eskandari–Murty–Nemoto construct a two-dimensional mixed motive over \(\mathbf Q\) whose period matrix contains Catalan's constant. Their work is explicitly designed around the philosophy that a sufficiently small geometric realization of a period may expose arithmetic information inaccessible in larger hypergeometric or modular realizations.

A basic period of their construction is

\[
\boxed{
G=
\iint_{\Delta}
\frac{dx\,dy}{1-x^2-y^2},
}
\]
where
\[
\Delta=\{(x,y):x\ge0,\ y\ge0,\ x+y\le1\}.
\]

This representation is particularly attractive because the associated period space is intrinsically two-dimensional: the geometry has already projected away the extraneous logarithmic and \(\pi\)-period directions that appear in more elementary hypergeometric descriptions.

The natural question is therefore:

> Can this two-dimensional Catalan motive be promoted to an Apéry or arithmetic-holonomy realization strong enough to prove \(\dim_{\mathbf Q}\langle1,G\rangle=2\)?

The calculations below investigate this question.

---

# 2. A moving Catalan period

Consider

\[
\boxed{
\Phi(z)
=
\iint_\Delta
\frac{dx\,dy}{1-z(x^2+y^2)}.
}
\]

Then

\[
\Phi(1)=G.
\]

Expanding at \(z=0\),

\[
\Phi(z)=\sum_{n\ge0}a_nz^n,
\qquad
a_n=\iint_\Delta(x^2+y^2)^n\,dx\,dy.
\]

Direct integration gives

\[
\boxed{
a_n=
\frac{1}{2^{n+1}(n+1)}
\sum_{k=0}^n
\frac{\binom nk}{2k+1}.
}
\]

The moments satisfy the exact recurrence

\[
\boxed{
(n+3)(2n+5)a_{n+2}
-(n+2)(3n+5)a_{n+1}
+(n+1)^2a_n=0.
}
\]

Thus the moving EMN period is unexpectedly low-rank.

If
\[
\theta=z\frac{d}{dz},
\]
then

\[
\boxed{
\left[
(\theta+1)(2\theta+1)
-z(\theta+1)(3\theta+2)
+z^2(\theta+1)^2
\right]\Phi(z)
=
\frac12.
}
\]

It is convenient to put

\[
H(z)=z\Phi(z).
\]

The equation then collapses to

\[
\boxed{
z(z-2)H''+(z-1)H'
=
\frac1{2(z-1)}.
}
\]

The associated homogeneous equation is

\[
z(z-2)Y''+(z-1)Y'=0,
\]

so the pure part is elementary. One may take

\[
1,
\qquad
S(z)=2\arcsin\sqrt{\frac z2}
\]

as homogeneous periods.

Hence the moving Catalan object is a single extension of an elementary rank-two local system.

---

## 3. Explicit uniformization

Set

\[
z=1-\cos S.
\]

Then

\[
S=2\arcsin\sqrt{\frac z2}.
\]

Writing

\[
\mathscr H(S)=H(1-\cos S),
\]

one obtains

\[
\boxed{
\mathscr H'(S)
=
\frac12
\log\tan\left(\frac\pi4+\frac S2\right).
}
\]

Integrating the elementary Fourier expansion gives

\[
\boxed{
\mathscr H(S)
=
\sum_{m\ge0}
\frac{(-1)^m}{(2m+1)^2}
\left(1-\cos((2m+1)S)\right).
}
\]

At

\[
S=\frac\pi2
\qquad(z=1)
\]

this gives

\[
\mathscr H(\pi/2)=G.
\]

This elementary trigonometric uniformization turns out to control most of the interesting extra structure.

---

# 4. A cubic fold and an exact distribution relation

The natural degree-three self-map of the pure coordinate is induced by the triple-angle identity.

Since

\[
R_3(z)=1-\cos(3S),
\]

with \(z=1-\cos S\), one finds

\[
\boxed{
R_3(z)=z(3-2z)^2.
}
\]

Thus

\[
S(R_3(z))=3S(z)
\]

on a compatible local branch.

Now use

\[
\prod_{j=0}^{2}
\tan\left(x+\frac{j\pi}{3}\right)
=
-\tan(3x).
\]

Applied to \(\mathscr H'\), this yields

\[
\mathscr H'(3S)
+
\sum_{j=0}^{2}
\mathscr H'
\left(S+\frac{2\pi j}{3}\right)
=0.
\]

After integration,

\[
\boxed{
\mathscr H(3S)
+
3\sum_{j=0}^{2}
\mathscr H
\left(S+\frac{2\pi j}{3}\right)
=
10G.
}
\]

The constant follows immediately from the Fourier series: the cubic trace isolates the odd integers divisible by \(3\), and

\[
\chi_{-4}(3)=-1.
\]

The factor

\[
1-\frac{\chi_{-4}(3)}{3^2}
=
1+\frac19
=
\frac{10}{9}
\]

is therefore built directly into the cubic geometry.

This gives a conceptual explanation for the ubiquitous conductor-\(12\) normalization \(10/9\).

---

# 5. The special point \(z=3/2\)

Take

\[
S=\frac{2\pi}{3}.
\]

Then

\[
z=1-\cos(2\pi/3)=\frac32.
\]

The cubic distribution gives

\[
\boxed{
\operatorname{PV}\Phi(3/2)
=
\frac{10}{9}G.
}
\]

Thus Catalan's constant admits a qualitatively different realization at a point outside the original Taylor disk:

\[
\boxed{
\operatorname{PV}
\iint_\Delta
\frac{dx\,dy}
{1-\frac32(x^2+y^2)}
=
\frac{10}{9}G.
}
\]

There is also the involution

\[
S\longmapsto\pi-S,
\]

which in the \(z\)-coordinate becomes

\[
z\longmapsto2-z.
\]

The Fourier expansion gives

\[
\boxed{
H(z)+H(2-z)=2G.
}
\]

Combining this with the value at \(3/2\),

\[
H(3/2)=\frac53G,
\]

gives

\[
H(1/2)=\frac13G,
\]

and therefore

\[
\boxed{
G=\frac32\Phi(1/2).
}
\]

This is notable because \(z=1/2\) lies strictly inside the Taylor disk.

Hence

\[
\boxed{
G=
\frac32
\sum_{n\ge0}\frac{a_n}{2^n}
}
\]

is an absolutely convergent rational Taylor realization derived directly from the EMN motive.

---

# 6. A quadratic-unit form

The same uniformization gives

\[
\varepsilon=2-\sqrt3,
\qquad
\varepsilon^{-1}=2+\sqrt3,
\]

and an exact representation

\[
\boxed{
G=
\frac32\,\operatorname{Ti}_2(2-\sqrt3)
+
\frac{\pi}{8}\log(2+\sqrt3),
}
\]

where

\[
\operatorname{Ti}_2(x)
=
\sum_{m\ge0}
\frac{(-1)^mx^{2m+1}}{(2m+1)^2}.
\]

The feature of interest here is that \(2-\sqrt3\) is a quadratic unit. Its small real embedding gives exponential decay without rational denominators; the conjugate embedding absorbs the corresponding growth.

This is a natural candidate for a number-field version of the arithmetic-holonomy programme, although the natural dilogarithmic companion still appears to require two LCM layers.

---

# 7. A \(3\)-adic Catalan realization

The same value \(z=3/2\) is especially interesting \(3\)-adically.

The Taylor series

\[
\Phi(3/2)
=
\sum_{n\ge0}a_n\left(\frac32\right)^n
\]

is divergent in the ordinary real Taylor sense, but the factor \(3^n\) causes it to converge \(3\)-adically.

Exact rational computation gives a \(3\)-adic number whose first thirty base-\(3\) digits agree with the conductor-\(12\) value

\[
L_3(2,\chi_{12}).
\]

Thus the calculations strongly identify

\[
\boxed{
\Phi_3(3/2)
=
L_3(2,\chi_{12}).
}
\]

This identity is also exactly what is predicted by the standard cubic distribution/Euler-factor relation for the \(3\)-adic dilogarithm. A fully formal proof should be written separately in the appropriate Coleman/syntomic normalization.

Combining the real and \(3\)-adic evaluations gives the adelic correspondence

\[
\boxed{
G
\quad\longleftrightarrow\quad
\frac9{10}L_3(2,\chi_{12}).
}
\]

This is significant because it supplies a genuinely non-Eisenstein realization of the same \(3\)-adic Catalan class appearing in the conductor-\(12\) modular construction.

---

# 8. Higher-pole EMN periods

The fixed EMN motive admits a large family of higher-pole periods.

For even \(m\ge t\ge0\), define

\[
I_{m,t}
=
\iint_\Delta
\frac{x^my^m}
{(1-x^2-y^2)^{t+1}}
\,dx\,dy.
\]

Write

\[
I_{m,t}=a_{m,t}+b_{m,t}G.
\]

For even \(t\ge2\), pole lowering gives

\[
\boxed{
I_{m,t}
=
-\frac{1}
{
2^{t+1}t(t-1)
\binom{2(m-t)}{m-t}
}
+
\frac{(m-1)^2}{4t(t-1)}
I_{m-2,t-2}.
}
\]

For example,

\[
I_{4,2}
=
-\frac{49}{384}
+\frac9{64}G.
\]

Iterating the recurrence gives a particularly simple formula for the Catalan coefficient:

\[
\boxed{
b_{m,t}
=
\frac1{4^m}
\binom m{m/2}
\binom mt.
}
\]

This closed form is one of the most useful outcomes of the calculation.

---

# 9. The central higher-pole family

Take

\[
m=4n,
\qquad
t=2n.
\]

Then

\[
I_n=
\iint_\Delta
\frac{x^{4n}y^{4n}}
{(1-x^2-y^2)^{2n+1}}
\,dx\,dy
=
a_n+b_nG
\]

with

\[
\boxed{
b_n=
\frac{\binom{4n}{2n}^2}{256^n}.
}
\]

Thus

\[
U_n:=256^nb_n=\binom{4n}{2n}^2\in\mathbf Z.
\]

Moreover

\[
b_n\sim\frac{1}{2\pi n}.
\]

The real saddle is

\[
\max_\Delta
\frac{x^4y^4}
{(1-x^2-y^2)^2}
=
\frac1{64},
\]

so

\[
I_n=64^{-n+o(n)}.
\]

Consequently

\[
-\frac{a_n}{b_n}
=
G+64^{-n+o(n)}.
\]

This is already an Apéry-shaped approximation.

Exact calculations further indicate that, after the normalization \(256^n\), the primary and rational companion belong to the same second-order holonomic system. Its asymptotic characteristic roots are

\[
\boxed{256,\qquad4.}
\]

Thus the dominant solution grows like \(256^n\), while the Catalan linear form grows like \(4^n\).

This explains both the quality and the limitation of the construction: after the natural integral normalization the recessive solution still grows rather than tends to zero.

---

# 10. Saddle killing

The natural idea is therefore to modify the numerator while retaining the homogeneous EMN reduction.

A first useful family is

\[
\boxed{
L_n(A,C)
=
\iint_\Delta
\frac{
x^{An}y^{An}(x^2-y^2)^{2Cn}
}
{(1-x^2-y^2)^{2n+1}}
\,dx\,dy.
}
\]

The factor

\[
(x^2-y^2)^{2Cn}
\]

kills the diagonal saddle \(x=y\) without introducing lower total degrees.

For \(D=1\), the exponential saddle can be evaluated exactly:

\[
\boxed{
M(A,C)
=
4^{-(A-1)}
\left(\frac{C}{A+C-2}\right)^C
\left(\frac{A-2}{A+C-2}\right)^{A-2}.
}
\]

Thus

\[
L_n(A,C)=M(A,C)^{\,n+o(n)}.
\]

The row

\[
(A,C)=(4,1)
\]

already gives

\[
M=\frac1{432}.
\]

More strikingly, the row

\[
\boxed{(A,C)=(10,3)}
\]

gives

\[
\boxed{
-\log M
=
18.9221280514\ldots.
}
\]

The corresponding \(G\)-coefficient has, in the computed range, only powers of \(2\) in its denominator. The observed primary dyadic cost is compatible with

\[
2^{26n+O(\log n)}b_n\in\mathbf Z.
\]

If this were the whole arithmetic cost, the margin would be

\[
18.9221280514-26\log2
=
\boxed{0.9003013569\ldots>0}.
\]

This was the first positive raw irrationality score found in the search.

However, the rational coefficient destroys the argument.

After clearing its odd denominator and then dividing by the full gcd of the two resulting integer coefficients, the primitive arithmetic cost is still vastly larger than the real decay—by more than \(20n\) in the tested range.

Thus the positive primary score does **not** extend to a positive primitive linear-form score.

---

# 11. Radial Padé factors

A different way to suppress the saddle is to use

\[
g=1-x^2-y^2
\]

and insert a polynomial in \(g^2\).

A particularly clean example is

\[
\boxed{
L_n=
\iint_\Delta
\frac{
x^{4n}y^{4n}(1-4g^2)^{2n}
}
{g^{2n+1}}
\,dx\,dy.
}
\]

The Catalan coefficient again satisfies

\[
\boxed{
256^nb_n\in\mathbf Z.
}
\]

For fixed \(g\),

\[
\max_{\Delta,\ 1-x^2-y^2=g}x^4y^4
=
\frac{\min(g,1-g)^4}{16},
\]

so introduce

\[
W(g)
=
\frac{\min(g,1-g)^4}{16g^2}.
\]

One obtains the exact weighted maximum

\[
\boxed{
\max_{0<g<1}
W(g)(1-4g^2)^2
=
\frac1{432}.
}
\]

Therefore

\[
\boxed{
0<256^nL_n
\le
G\left(\frac{16}{27}\right)^n.
}
\]

This is a genuine “free integration” phenomenon:

\[
\frac{16}{27}<1.
\]

After normalizing the \(G\)-coefficient integrally, the real linear form already decreases geometrically.

Unfortunately the rational companion again acquires large odd denominators. Increasing the radial degree improves the real weighted norm but systematically makes this arithmetic defect worse.

---

# 12. Jacobi and Racah boundary cancellation

The EMN pole-lowering operation produces explicit rational boundary terms. This makes a Padé strategy possible: choose the numerator so that these boundary functionals vanish.

A simple finite-difference construction is

\[
\boxed{
F_n=
x^{2n}y^{2n}
\left((1-x^2)^n+(1-y^2)^n\right).
}
\]

When expanded, the relevant coefficients are

\[
c_0=2,
\qquad
c_k=(-1)^k\binom nk.
\]

The \(s\)-th boundary functional is a polynomial of degree \(s<n\) in \(k\), hence

\[
\sum_{k=0}^n(-1)^k\binom nkP_s(k)=0.
\]

Therefore this numerator annihilates all \(n\) boundary terms exactly.

A one-parameter version is

\[
\boxed{
\begin{aligned}
F_{n,\beta}
={}&x^{2n}y^{2n}\Big[
(1-x^2)^{n-1}
\{\beta-(\beta+n)x^2\}\\
&\qquad+
(1-y^2)^{n-1}
\{\beta-(\beta+n)y^2\}
\Big],
\end{aligned}
}
\]

which kills the first \(n-1\) boundary terms.

The terminal numerator is hypergeometric/Jacobi:

\[
T^n[x^{2n}(1-x^2)^n]
\propto
{}_2F_1
\left(-n,n+\frac12;\frac12;x^2\right),
\]

equivalently a Jacobi polynomial

\[
P_n^{(-1/2,0)}(1-2x^2).
\]

A balanced constant-degree construction instead leads to quadratic-lattice orthogonality and Racah/dual-Hahn coefficients.

For instance, on a suitable shell the boundary annihilator is

\[
\boxed{
c_{n,k}
=
(-1)^k
\binom{2n}{k}
\binom{4n}{2k}.
}
\]

Thus the Padé structure of the EMN motive naturally produces classical discrete orthogonal polynomials.

This is conceptually attractive, but the numerical outcome is again negative: the Jacobi family has cheap coefficients but insufficient real decay; the balanced Racah family has strong Catalan cancellation but exponential coefficient height which is too large.

---

# 13. Cubic conditional holonomy

The cubic map also produces a conditional function directly.

Since

\[
R_3(1)=1,
\qquad
R_3'(1)=-3,
\]

the logarithmic monodromy at \(z=1\) cancels in

\[
H(R_3(z))+3H(z).
\]

Its value at \(z=1\) is \(4G\). Hence, if \(G\in\mathbf Q\),

\[
\boxed{
K(z)=H(R_3(z))+3H(z)-4G
}
\]

has rational Taylor coefficients and is holomorphic across the Catalan fold \(z=1\).

Its coefficients exhibit a striking dyadic pattern:

\[
v_2([z^n]K)=n+O(\log n)
\]

in the computed range.

Thus the conditional generator has a genuine \(2\)-adic slope

\[
\boxed{\varsigma_2=1.}
\]

Unfortunately \(R_3^{-1}(1)\) contains

\[
1,\qquad
1-\frac{\sqrt3}{2},\qquad
1+\frac{\sqrt3}{2},
\]

so cancelling only the branch at \(1\) introduces a closer singularity at

\[
1-\frac{\sqrt3}{2}.
\]

Trying to cancel all three simultaneously gives a rigidity phenomenon: the unique local coefficient pattern is exactly the full cubic trace

\[
(1,3,3,3),
\]

which collapses to the constant \(10G\).

Thus there is no nonconstant full-orientation cubic conditional generator.

---

# 14. Negative conclusions

Several attractive routes can now be ruled out, at least within the classes tested.

### 14.1 The naive moving period

The function \(\Phi(z)-G\) does not become holomorphic at \(z=1\) under the assumption \(G\in\mathbf Q\): the logarithmic coefficient is independent of \(G\).

### 14.2 Parameter-derivative Gauss systems

Quarter-parameter zero-balanced Gauss deformations can isolate \(G\) in a connection coefficient, but the unwanted \(\log2\) or second-LCM contribution reappears elsewhere.

### 14.3 Pure radial saddle engineering

Real decay can be made very strong, but extra radial degree increases the rational companion denominator too rapidly.

### 14.4 Jacobi/Racah boundary cancellation

The boundary terms can be cancelled exactly, but the remaining real/arithmetic score remains negative.

### 14.5 Homogeneous saddle killers

These can produce positive primary margins, but their primitive rational companions remain too expensive.

### 14.6 Odd-pole reduction

A census using the distinct odd-pole operator \(D_1+D_2\) also gives negative primitive scores.

### 14.7 Full cubic trace

Complete cancellation of cubic-fold singularities is rigid and gives only the constant distribution relation.

---

# 15. What appears genuinely new

The most useful discoveries are therefore structural rather than an irrationality theorem.

### A. A very small moving Picard–Fuchs model

\[
\boxed{
z(z-2)H''+(z-1)H'=\frac1{2(z-1)}
}
\]

gives a rank-three extension model with elementary rank-two pure part.

### B. Exact cubic Catalan distribution

\[
\boxed{
\mathscr H(3S)
+
3\sum_{j=0}^{2}
\mathscr H(S+2\pi j/3)
=
10G.
}
\]

This geometrically explains the factor \(10/9\).

### C. A non-Eisenstein \(3\)-adic Catalan realization

The EMN period at \(z=3/2\) realizes, at the two places,

\[
\boxed{
\frac{10}{9}G
\quad\leftrightarrow\quad
L_3(2,\chi_{12}),
}
\]

at least with overwhelming computational evidence and the expected standard \(p\)-adic distribution law.

### D. Explicit higher-pole Catalan coefficient

\[
\boxed{
b_{m,t}
=
4^{-m}
\binom m{m/2}\binom mt.
}
\]

### E. Genuine EMN Apéry systems

The central higher-pole family gives a dominant/recessive structure with characteristic rates

\[
256,\qquad4.
\]

### F. A true primary free-integration family

\[
256^n
\iint_\Delta
\frac{x^{4n}y^{4n}(1-4g^2)^{2n}}
{g^{2n+1}}dx\,dy
=
O\left((16/27)^n\right).
\]

### G. Orthogonal-polynomial structure

The natural Padé problem on the EMN motive produces Jacobi and Racah/dual-Hahn systems rather than arbitrary coefficient arrays.

These observations suggest that the Catalan motive is considerably more structured than its classical rational approximations alone reveal.

---

# 16. The remaining obstruction

Every promising construction eventually encounters the same arithmetic phenomenon.

At the level of the Catalan coefficient \(b_n\), substantial cancellation occurs:

- often only powers of \(2\) remain;
- sometimes the normalized real linear form already decays geometrically;
- the motive itself has only dimension \(2\).

But the rational coefficient \(a_n\) retains a much larger denominator lattice.

Thus the current obstruction can be summarized as

\[
\boxed{
\text{Catalan extension class has good arithmetic;}
\quad
\text{its chosen rational splitting does not.}
}
\]

This suggests that further progress should focus not on improving the real saddle but on changing the **integral realization of the rational splitting**.

That is precisely the philosophy behind arithmetic holonomy, syntomic realizations and multi-lattice Apéry theory.

---

# 17. Proposed next programme

The work points to four concrete directions.

### 17.1 Compute the integral structure of the EMN motive itself

Rather than bounding denominators after writing
\[
I_n=a_n+b_nG,
\]
one should determine the natural integral Betti/de Rham lattice of the two-dimensional motive and express the higher-pole periods in that lattice.

The enormous denominators in \(a_n\) may partly be artefacts of the coordinate \(1,G\), rather than intrinsic indices of the motive.

### 17.2 Construct the syntomic/Frobenius realization

The identity at \(z=3/2\) suggests that the EMN motive carries the same \(3\)-adic Catalan extension as the conductor-\(12\) modular object.

An explicit Frobenius matrix for the EMN rank-two extension would allow the real and \(3\)-adic realizations to be compared directly, without building a separate rational approximating sequence.

### 17.3 Search for an adelic projector

The ideal construction would be an algebraic projector
\[
P
\]
which simultaneously:

- retains the Catalan extension;
- kills the expensive rational splitting;
- has small archimedean norm;
- acquires positive \(2\)- or \(3\)-adic Frobenius slope.

The Jacobi/Racah calculations show that the necessary projectors already have a classical orthogonal-polynomial structure.

### 17.4 Compare the EMN and conductor-\(12\) motives

The numerical and distribution-theoretic identity

\[
G
\leftrightarrow
\frac9{10}L_3(2,\chi_{12})
\]

should be upgraded to an explicit morphism or comparison of extension classes.

If such a comparison exists motivically or syntomically, it would give a precise example of the general “one class, two worlds” principle:

\[
\boxed{
\text{one extension class}
\quad
\rightsquigarrow
\quad
\text{different integral realizations}
\quad
\rightsquigarrow
\quad
\text{different arithmetic slopes}.
}
\]

That may ultimately be more useful than any individual family of rational approximants.

---

# 18. Conclusion

The Eskandari–Murty–Nemoto motive provides the qualitatively different Catalan realization that the earlier modular programme was missing. Their 2025 preprint explicitly constructs a two-dimensional mixed motive over \(\mathbf Q\) having \(G\) as a period and derives linear forms in \(1,G\).

Developing that realization produced:

\[
\boxed{
\begin{gathered}
\text{an elementary Picard--Fuchs extension},\\
\text{a cubic distribution law},\\
\text{a conductor-\(12\) real/\(3\)-adic bridge},\\
\text{new Apéry-type higher-pole systems},\\
\text{Jacobi/Racah Padé structures},\\
\text{and explicit free-integration phenomena}.
\end{gathered}}
\]

None yet proves Catalan's irrationality.

The sharp remaining problem is no longer to find a rapidly decaying Catalan period. Such periods now exist. It is to find an **integral realization in which the rational companion loses approximately one denominator layer**.

In that form, the problem is highly compatible with the wider arithmetic-holonomy programme: the geometry, real decay, \(p\)-adic realization and low motivic rank are already present. The missing ingredient is the correct lattice.

## Reference

P. Eskandari, V. Kumar Murty and Y. Nemoto, *Mixed motives and linear forms in the Catalan constant*, arXiv:2510.20648, 2025. The authors construct a two-dimensional mixed motive over \(\mathbf Q\) having Catalan's constant as a period and derive explicit families of linear forms in \(1\) and \(G\).