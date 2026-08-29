# A New Mixed-Motive and Cyclotomic Apéry Programme for Catalan’s Constant

## Research report

### Abstract

This report records an attempt to find a qualitatively different realization of Catalan's constant
\[
G=\beta(2)=\sum_{n\ge0}\frac{(-1)^n}{(2n+1)^2}
\]
that might evade the denominator and holonomy obstructions encountered by the modular/Eisenstein Catalan constructions developed previously.

Three related but distinct realizations were investigated:

1. a signature-four hypergeometric/symmetric-square realization;
2. direct cyclotomic two-dimensional Beukers integrals;
3. the new two-dimensional mixed motive of Eskandari–Murty–Nemoto (EMN), together with moving, higher-pole and Padé deformations of that motive.

The main conclusions are as follows.

First, the signature-four construction genuinely changes the Hodge placement of \(G\) and produces a one-LCM parameter-derivative extension, but its local connection matrix contains unavoidable nuisance periods and does not give a rationality-dependent holonomy fold.

Second, a direct cyclotomic cellular integral produces a new-looking second-order Apéry recurrence converging to \(G\). This gives a useful new Catalan Apéry system, but its recessive root and companion denominators remain too large for irrationality.

Third, the EMN motive is the most promising new geometric host. Its higher-pole reduction can be made completely explicit. A radial Padé modification produces a genuine "free integration" phenomenon: after the natural integral normalization, the real linear form already decays geometrically. The remaining obstruction is entirely arithmetic, namely the denominator of the rational companion.

Several exact Jacobi/Racah Padé constructions were found that annihilate the rational boundary contributions created during pole reduction. None yet combines sufficient real decay with sufficiently cheap arithmetic. The remaining search has therefore become a finite and explicit optimization problem rather than an open-ended hunt for new representations of \(G\).

A further adelic test shows that the natural EMN approximants possess a strong \(2\)-adic slope but apparently no \(3\)-adic slope. The highest-priority next experiment is to test the new boundary-cancelled and radial-Padé EMN realizations for a positive \(3\)-adic slope.

---

# 1. Motivation

Previous work on Catalan's constant had reached a fairly robust obstruction.

The standard modular realization places \(G\) in a weight-three Eisenstein/Eichler-integral architecture. Its arithmetic naturally carries approximately two LCM layers. Attempts involving Atkin–Lehner transformations, \(\mu_4\)-polylogarithmic hosts, additional period classes and ordinary arithmetic holonomy did not overcome that denominator cost.

The archive summarized the remaining qualitative escape route as:

\[
\boxed{\text{find a non-Eisenstein realization of the same Catalan class}.}
\]

The present work pursued that possibility.

The guiding principle throughout was that

\[
\text{same period}
\neq
\text{same motive}
\neq
\text{same integral realization}
\neq
\text{same Diophantine cost}.
\]

Thus merely finding another identity for \(G\) is of little value. The new realization must alter at least one of:

- the Hodge/integration depth;
- the integral lattice;
- the singularity geometry;
- the archimedean decay;
- the \(p\)-adic Frobenius slopes.

---

# 2. Signature-four elliptic realization

Let
\[
F(t)={}_2F_1\!\left(\frac14,\frac34;1;t\right).
\]

A first qualitatively different realization is

\[
\boxed{
\int_0^1F(t)^2\,dt=\frac{16G}{\pi^2}.
}
\]

Thus \(G\) occurs in the symmetric square of a rank-two signature-four elliptic local system.

Define
\[
I_n=\int_0^1t^nF(t)^2\,dt.
\]

The symmetric-square differential equation is

\[
\boxed{
8t^2(1-t)^2Y'''
+24t(1-t)(1-2t)Y''
+2(27t^2-27t+4)Y'
+3(2t-1)Y=0.
}
\]

Integrating against \(t^n\) yields

\[
\boxed{
\begin{aligned}
8n^3I_{n-1}
&-(2n+1)(8n^2+8n+3)I_n\\
&+2(n+1)(2n+1)(2n+3)I_{n+1}
=\frac8{\pi^2}.
\end{aligned}
}
\]

In particular,

\[
I_0=\frac{16G}{\pi^2},
\qquad
I_1=\frac{8G+\frac43}{\pi^2}.
\]

Hence

\[
\boxed{
\frac1{\pi^2}
=
\frac34\left(I_1-\frac12I_0\right).
}
\]

Thus both \(G/\pi^2\) and \(1/\pi^2\) belong to the same rank-three period system.

This is considerably better structurally than an isolated formula for \(G/\pi^2\).

## Status

The differential equation, moment recurrence and first moment evaluations were derived exactly.

---

# 3. Parameter derivative and a one-LCM extension

Consider the zero-balanced family

\[
F_a(t)={}_2F_1(a,1-a;1;t)
\]

and its parameter derivative

\[
D(t)=
\left.
\frac{\partial}{\partial a}F_a(t)
\right|_{a=1/4}.
\]

Differentiating the Gauss equation gives

\[
\boxed{
t(1-t)D''+(1-2t)D'-\frac3{16}D
=\frac12F.
}
\]

Thus \(D\) is an extension of the signature-four elliptic system.

After the natural global-boundedness scaling,

\[
F(64x)
=
\sum_{n\ge0}
\binom{4n}{2n}\binom{2n}{n}x^n.
\]

The derivative coefficients are

\[
[x^n]D(64x)
=
\binom{4n}{2n}\binom{2n}{n}
\left[
4\sum_{k=0}^{n-1}
\left(
\frac1{4k+1}-\frac1{4k+3}
\right)
\right].
\]

The notable arithmetic feature is

\[
\boxed{
d_n[x^n]D(64x)\in\mathbf Z,
\qquad
d_n=\operatorname{lcm}(1,\ldots,n).
}
\]

Thus the extension itself has only one LCM layer.

The connection formula contains

\[
\psi_1\!\left(\frac14\right)
-
\psi_1\!\left(\frac34\right)
=16G.
\]

Consequently Catalan's constant appears as a connection coefficient of this arithmetic differential system.

This was the first indication that the usual two-integration Catalan cost is not intrinsic to the period \(G\).

## Obstruction

Unfortunately, \(G\) belongs to the regular part of the connection matrix rather than to the first logarithmic residue. Rationality of \(G\) therefore does not remove the local logarithmic monodromy.

Attempts to close the corresponding signature-four quadratic isogeny on the one- or two-dimensional normal-function space also fail: additional de Rham directions appear.

The construction is therefore valuable as an explanatory model but does not by itself supply the required conditional holonomy fold.

---

# 4. A direct cyclotomic cellular realization

A genuinely different elementary period is

\[
\boxed{
G=
\int_0^1\!\!\int_0^1
\frac{dx\,dy}{1+x^2y^2}.
}
\]

For integrals of the form

\[
\int_0^1\!\!\int_0^1
\frac{x^Ay^B}{(1+x^2y^2)^r}\,dx\,dy,
\]

the substitution \(u=xy\) reduces everything to the one-variable families

\[
J_m^{(r)}
=
\int_0^1
\frac{u^m}{(1+u^2)^r}\,du,
\]

and

\[
K_m^{(r)}
=
-\int_0^1
\frac{u^m\log u}{(1+u^2)^r}\,du.
\]

These satisfy elementary rational recurrences.

This gives an exact period decomposition in the basis

\[
1,\quad \pi,\quad\log2,\quad G,\quad\pi^2.
\]

Congruence conditions on the monomial exponents can force some or all of the unwanted period coordinates to vanish.

---

# 5. A new cyclotomic Apéry recurrence

The particularly clean family

\[
\boxed{
C_n=
\int_0^1\!\!\int_0^1
\frac{
x^{2n}(1-x^4)^n
y^{2n}(1-y^4)^n
}
{(1+x^2y^2)^{2n+1}}
\,dx\,dy
}
\]

has period space exactly

\[
C_n=A_nG+B_n,
\qquad
A_n,B_n\in\mathbf Q.
\]

All \(\pi\), \(\log2\), and \(\pi^2\) coordinates cancel structurally.

Define

\[
U_n=(-1)^n16^nA_n,
\qquad
V_n=(-1)^n16^nB_n.
\]

Experimentally,

\[
U_0=1,\quad
U_1=28,\quad
U_2=2596,\quad
U_3=311536,\quad
U_4=41759524,\ldots
\]

and every computed \(U_n\) is integral.

Both \(U_n\) and \(V_n\) satisfy the same recurrence

\[
\boxed{
\begin{aligned}
&(n+2)^2(2n+3)^2(20n^2+32n+13)X_{n+2}\\
&\quad
-4P(n)X_{n+1}\\
&\quad
-256(n+1)^2(2n+1)^2
(20n^2+72n+65)X_n=0,
\end{aligned}}
\]

where

\[
\begin{aligned}
P(n)=
&3520n^6+26752n^5+83024n^4\\
&+134592n^3+120196n^2+56088n+10699.
\end{aligned}
\]

The initial companion data are

\[
V_0=0,\qquad V_1=-26.
\]

Consequently

\[
-\frac{V_n}{U_n}\longrightarrow G.
\]

## Status

This recurrence was discovered by exact guessing and verified on 61 exact terms.

A creative-telescoping or Picard–Fuchs proof remains to be supplied before treating it as a theorem.

---

# 6. Why this new recurrence does not prove irrationality

The characteristic equation at infinity is

\[
\lambda^2-176\lambda-256=0,
\]

with roots

\[
\boxed{
\lambda_\pm=88\pm40\sqrt5.
}
\]

The recessive root is

\[
\lambda_-=88-40\sqrt5
=-1.442719\ldots.
\]

Hence

\[
|\lambda_-|>1.
\]

So even before clearing the remaining companion denominators, the naturally normalized recessive solution does not tend to zero.

Furthermore, experimentally,

\[
\boxed{
\operatorname{den}(V_n)\mid d_{2n}^{\,2}
}
\]

through the exact range tested.

Thus this row reproduces the familiar weight-two denominator obstruction in a new geometric guise.

The importance of the recurrence is therefore conceptual rather than immediately Diophantine: it shows that a very simple cyclotomic cell generates a genuine Apéry system for \(G\).

---

# 7. The Eskandari–Murty–Nemoto mixed motive

The decisive external development is the recent construction of Eskandari, Murty and Nemoto.

They construct a two-dimensional mixed motive over \(\mathbf Q\) having \(G\) as a period and produce a large family of linear forms in \(1\) and \(G\).

Their basic geometric period is

\[
\boxed{
G=
\iint_\Delta
\frac{dx\,dy}{1-x^2-y^2},
}
\]

where

\[
\Delta
=
\{(x,y):x\ge0,\ y\ge0,\ x+y\le1\}.
\]

The order-four symmetry

\[
\sigma(x,y)=(-y,x)
\]

is used to project the ambient relative cohomology onto a two-dimensional Catalan motive.

In an appropriate basis the period matrix has the form

\[
\boxed{
\begin{pmatrix}
1&G\\
0&-\pi(2\pi i)
\end{pmatrix}.
}
\]

Over \(\mathbf Q(i)\) the graded pieces become Tate.

This is exactly the kind of non-Eisenstein realization that the previous Catalan obstruction analysis had identified as missing.

---

# 8. The naive moving EMN period

The most obvious deformation is

\[
\Phi(z)
=
\iint_\Delta
\frac{dx\,dy}
{1-z(x^2+y^2)},
\qquad
\Phi(1)=G.
\]

Write

\[
\Phi(z)=\sum_{n\ge0}a_nz^n.
\]

The moments are

\[
\boxed{
a_n=
\frac1{2^{n+1}(n+1)}
\sum_{k=0}^n
\frac{\binom nk}{2k+1}.
}
\]

They satisfy

\[
\boxed{
(n+3)(2n+5)a_{n+2}
-(n+2)(3n+5)a_{n+1}
+(n+1)^2a_n=0.
}
\]

Putting

\[
H(z)=z\Phi(z)
\]

gives the remarkably small Picard–Fuchs equation

\[
\boxed{
z(z-2)H''+(z-1)H'
=
\frac1{2(z-1)}.
}
\]

The homogeneous equation

\[
z(z-2)H''+(z-1)H'=0
\]

has solutions

\[
1,
\qquad
2\arcsin\sqrt{\frac z2}.
\]

Moreover,

\[
\boxed{
H'(z)
=
\frac{
\operatorname{artanh}\sqrt{z/(2-z)}
}{
\sqrt{z(2-z)}
}.
}
\]

This gives an explicit rank-three Gauss–Manin system.

---

# 9. Structural obstruction to the naive EMN deformation

Near \(z=1\),

\[
\boxed{
H(z)
=
G+\frac{1-z}{2}\log(1-z)+O(1-z).
}
\]

The logarithmic branch coefficient is independent of \(G\).

Therefore the hypothesis

\[
G\in\mathbf Q
\]

does not remove the branch.

This kills the simplest conditional-holonomy strategy.

There is also a quantitative version.

After the natural scaling \(z=2x\), the coefficients have denominator type approximately

\[
e=1,\qquad b=2.
\]

The closest logarithmic singularity occurs at

\[
x=\frac12.
\]

Any simply connected continuation disk consequently has conformal radius at most \(2\), while the elementary flat arithmetic-holonomy cost is already at least

\[
\tau^\flat=\frac32.
\]

Thus

\[
\log2<\frac32.
\]

The Bost–Charles/CDT entry inequality cannot be reached, even before the positive integration correction is added.

Hence

\[
\boxed{
\text{the naive moving EMN period is not a viable CDT host.}
}
\]

This is a genuine structural no-go result for this specific deformation.

---

# 10. Exact higher-pole reduction

The fixed EMN motive remains much more promising.

For even \(m\ge t\ge0\), define

\[
I_{m,t}
=
\iint_\Delta
\frac{x^my^m}
{(1-x^2-y^2)^{t+1}}\,dx\,dy.
\]

Then

\[
I_{m,t}=a_{m,t}+b_{m,t}G.
\]

For \(t\ge2\), direct pole lowering gives

\[
\boxed{
I_{m,t}
=
-\frac{
1
}{
2^{t+1}t(t-1)
\binom{2(m-t)}{m-t}
}
+
\frac{(m-1)^2}{4t(t-1)}
I_{m-2,t-2}.
}
\]

As a check,

\[
\boxed{
I_{4,2}
=
-\frac{49}{384}+\frac9{64}G,
}
\]

agreeing with the EMN example.

The Catalan coefficient has a particularly simple closed form:

\[
\boxed{
b_{m,t}
=
\frac1{4^m}
\binom m{m/2}\binom mt.
}
\]

This formula is one of the most useful results of the calculation.

---

# 11. The central-pole regime

Take

\[
m=4n,
\qquad
t=2n.
\]

Then

\[
\boxed{
b_n
=
\frac{\binom{4n}{2n}^{\,2}}{4^{4n}}
\sim\frac1{2\pi n}.
}
\]

Thus the coefficient of \(G\) is only polynomially small.

Meanwhile the positive integral has exponential kernel

\[
\frac{x^4y^4}{(1-x^2-y^2)^2}.
\]

Its maximum on the triangle occurs at

\[
x=y=\frac12
\]

and equals

\[
\frac1{64}.
\]

Therefore

\[
\boxed{
I_{4n,2n}
=
64^{-n+o(n)}.
}
\]

This is exactly the desired qualitative shape:

\[
\frac{-a_n}{b_n}
=
G+
64^{-n+o(n)}.
\]

The obstruction is purely arithmetic: the rational coefficient \(a_n\) has an exponentially expensive denominator.

---

# 12. Radial Padé deformation

The central saddle occurs on the level

\[
g=1-x^2-y^2=\frac12.
\]

This suggests modifying the numerator by a polynomial in \(g^2\).

A particularly useful family is

\[
\boxed{
L_n=
\iint_\Delta
\frac{
x^{4n}y^{4n}
(1-4g^2)^{2n}
}
{g^{2n+1}}
\,dx\,dy
=
a_n+b_nG.
}
\]

The \(G\)-coefficient has excellent arithmetic:

\[
\boxed{
256^n b_n\in\mathbf Z.
}
\]

More explicitly,

\[
256^n b_n
=
\binom{4n}{2n}
\sum_{j=0}^{n}
(-4)^j
\binom{2n}{j}
\binom{4n}{2n-2j}.
\]

The restriction to \(j\le n\) reflects the fact that higher powers cancel the pole completely and contribute no \(G\)-term.

---

# 13. A genuine archimedean threshold crossing

For fixed

\[
g=1-x^2-y^2,
\]

one has

\[
\max_\Delta x^4y^4
=
\frac{\min(g,1-g)^4}{16}.
\]

Thus the one-dimensional weight is

\[
\boxed{
W(g)=
\frac{\min(g,1-g)^4}{16g^2}.
}
\]

For the radial factor,

\[
\boxed{
256
\max_{0<g<1}
W(g)(1-4g^2)^2
=
\frac{16}{27}.
}
\]

Consequently

\[
\boxed{
0<256^nL_n
\ll
\left(\frac{16}{27}\right)^n.
}
\]

This is significant.

After the normalization that makes the Catalan coefficient integral, the real linear form already tends geometrically to zero.

Thus

\[
\boxed{
\text{the archimedean side of an Apéry argument has been solved for this realization.}
}
\]

This is the clearest instance found here of the anticipated "free integration" phenomenon.

---

# 14. The remaining denominator wall

Unfortunately the rational companion does not share the excellent arithmetic of the \(G\)-coefficient.

Exact computations indicate a denominator law of approximately

\[
\operatorname{den}(256^na_n)
\mid
L_{12n}^{\,2}
\]

or a comparably expensive two-LCM structure.

The precise sharp theorem has not yet been proved; this statement should be regarded as an exact-data conjecture.

Its exponential cost is far larger than the \(16/27\) real gain.

Thus the radial construction demonstrates that good real geometry is possible, but it does not yet prove irrationality.

---

# 15. Degree-two radial optimization

The one-dimensional reformulation makes radial optimization very efficient.

For a linear real factor

\[
1-cg^2,
\]

the real optimum occurs near

\[
c\approx4.6240,
\]

for which numerically

\[
256\,
\|W(g)(1-cg^2)\|_\infty
\approx0.865<1.
\]

No integer \(c\) gives the same improvement: \(c=4\) lies at the critical boundary.

Searching integer quadratic factors produced the candidate

\[
\boxed{
q(h)=1-10h+23h^2,
}
\]

with numerical score

\[
\boxed{
256
\max W(g)|q(g^2)|
\approx0.50119.
}
\]

This gives an even larger real margin while keeping integral polynomial coefficients.

However, increasing radial degree also increases the rational companion denominator. In the families tested, that arithmetic growth overwhelms the real improvement.

Thus radial degree by itself is not the missing ingredient.

---

# 16. Pole-lowering Padé conditions

The rational denominator enters through explicit boundary terms during pole reduction.

This suggests cancelling those terms directly.

Let

\[
D=D_2D_1
\]

denote the relevant reduction operator.

The ideal Padé conditions are

\[
\boxed{
R_0(F_n)=
R_1(DF_n)=
\cdots=
R_{n-1}(D^{n-1}F_n)=0.
}
\]

If these hold, the successive rational boundary contributions disappear.

The remaining simple-pole numerator can additionally be constrained to have zero rational period.

This turns the numerator search into exact finite-dimensional linear algebra.

---

# 17. Jacobi finite-difference family

A particularly inexpensive basis gives an explicit solution of the boundary conditions.

Set

\[
\boxed{
F_n=
x^{2n}y^{2n}
\left[
(1-x^2)^n+(1-y^2)^n
\right].
}
\]

Expanding,

\[
c_0=2,
\qquad
c_k=(-1)^k\binom nk.
\]

After the \(s\)-th reduction, the corresponding boundary functional is a polynomial in \(k\) of degree \(s\).

Therefore for \(s<n\),

\[
\sum_{k=0}^n
(-1)^k\binom nkP_s(k)=0.
\]

Hence

\[
\boxed{
\text{all \(n\) pole-lowering boundary contributions vanish exactly.}
}
\]

A one-parameter version is

\[
\boxed{
\begin{aligned}
F_{n,\beta}
={}&x^{2n}y^{2n}
\Big[
(1-x^2)^{n-1}
\{\beta-(\beta+n)x^2\}\\
&+
(1-y^2)^{n-1}
\{\beta-(\beta+n)y^2\}
\Big].
\end{aligned}}
\]

This kills the first \(n-1\) boundary terms and leaves only the final one.

The construction is naturally Jacobi/hypergeometric rather than accidental finite-dimensional cancellation.

---

# 18. Racah/dual-Hahn structure

A balanced constant-degree shell produces another exact Padé system:

\[
\phi_k=
x^{2n+2k}y^{6n-2k}
+
x^{6n-2k}y^{2n+2k}.
\]

The pole-lowering boundary rows become polynomials in the quadratic lattice variable

\[
(n-k)^2.
\]

The unique annihilating weights are

\[
\boxed{
c_{n,k}
=
(-1)^k
\binom{2n}{k}
\binom{4n}{2k}.
}
\]

This is the discrete orthogonality structure characteristic of dual-Hahn/Racah systems.

Thus the EMN Padé problem naturally generates classical discrete orthogonal polynomials.

This is probably conceptually significant even independently of irrationality.

---

# 19. Why the first Jacobi/Racah constructions fail

The Jacobi family has small coefficients and excellent exact boundary cancellation, but its real saddle is too large.

The balanced Racah family gives much stronger Catalan cancellation, but its coefficient height becomes enormous.

In both cases the primitive integer height grows faster than the resulting approximation error shrinks.

Thus:

\[
\boxed{
\text{boundary cancellation alone is insufficient.}
}
\]

The remaining numerator should probably combine:

- a low-degree radial factor for real saddle suppression;
- a separated/Jacobi component for divisibility under pole lowering.

This motivates a hybrid construction rather than either mechanism alone.

---

# 20. A pure-\(G\) terminal condition

At simple pole, the EMN period map is

\[
H\longmapsto a_H+b_HG.
\]

Therefore a particularly attractive terminal Padé condition is

\[
\boxed{
a_H=0,\qquad b_H\ne0.
}
\]

Then the terminal simple-pole contribution is a pure multiple of \(G\).

For example,

\[
\boxed{
H=
1+569x^2y^2-2800x^4y^4
}
\]

satisfies

\[
\boxed{
\iint_\Delta
\frac{H\,dx\,dy}
{1-x^2-y^2}
=
\frac{13}{2}G.
}
\]

Thus pure-\(G\) polynomial sections exist explicitly.

The corresponding backwards-Padé problem is to pull such a numerator through \(D_2D_1\) while simultaneously controlling the earlier boundary terms.

This remains a promising direction.

---

# 21. Motivic splitting does not immediately prove irrationality

One tempting possibility was:

\[
G\in\mathbf Q
\quad\Longrightarrow\quad
\text{EMN extension splits}.
\]

This implication is false.

Over \(\mathbf Q(i)\), the motive is an extension involving Tate pieces, and the extension class is measured relative to the relevant Tate periods and Hodge filtration.

Rationality of the real off-diagonal entry \(G\) does not imply motivic splitting.

Thus non-splitting of the EMN motive alone is not an irrationality proof.

---

# 22. Adelic test

The previous Catalan programme identified a missing \(3\)-adic realization as potentially decisive.

A natural question is therefore whether the EMN rational approximants have useful \(p\)-adic convergence.

For

\[
I_n=a_n+b_nG
\]

define

\[
r_n=-\frac{a_n}{b_n}.
\]

The exact rational quantities

\[
r_n-r_{n-1}
\]

were examined \(p\)-adically.

For the natural central EMN family, the experimental result is:

\[
\boxed{
v_2(r_n-r_{n-1})
\sim 5.9n,
}
\]

while

\[
\boxed{
v_3(r_n-r_{n-1})
=O(1)
}
\]

over the range tested, through approximately \(n=60\).

Thus the natural EMN row appears to have

\[
\sigma_2\approx6,
\qquad
\sigma_3=0.
\]

## Status

These are computational observations, not proved asymptotic theorems.

The absence of a visible \(3\)-slope is nevertheless very clear numerically.

---

# 23. Overall assessment

The work has ruled out several superficially promising ideas:

- the naive signature-four isogeny fold;
- the naive moving EMN period;
- simple radial saddle suppression;
- pure Jacobi boundary cancellation;
- pure Racah boundary cancellation;
- motivic non-splitting as an immediate irrationality argument.

But it has also uncovered substantially new positive structure.

Most importantly:

\[
\boxed{
\text{there exist EMN realizations for which the naturally normalized real form already decays geometrically.}
}
\]

The example

\[
x^{4n}y^{4n}(1-4g^2)^{2n}
\]

gives the explicit normalized bound

\[
\left(\frac{16}{27}\right)^n.
\]

This means the old Catalan difficulty should no longer be described simply as "not enough analytic decay."

For this realization,

\[
\boxed{
\text{the only serious obstruction is the rational companion lattice.}
}
\]

That is a much narrower problem.

---

# 24. Current research target

The next numerator should simultaneously possess:

\[
\boxed{
\begin{array}{ll}
\textbf{A.}&
256^n\text{-integral or comparably cheap \(G\)-coefficient},\\[1mm]
\textbf{B.}&
\text{weighted real norm strictly below the irrationality threshold},\\[1mm]
\textbf{C.}&
\text{large content under repeated }D_2D_1,\\[1mm]
\textbf{D.}&
\text{cancellation of most rational boundary functionals},\\[1mm]
\textbf{E.}&
\text{preferably a positive \(3\)-adic slope}.
\end{array}
}
\]

The evidence suggests that no single simple polynomial mechanism achieves all five.

The most plausible architecture is hybrid:

\[
\boxed{
\text{radial Chebyshev factor}
\times
\text{Jacobi/Racah arithmetic factor}.
}
\]

The radial factor attacks the real saddle.

The Jacobi/Racah factor attacks the rational companion denominator.

These should be optimized jointly rather than sequentially.

---

# 25. Highest-priority next calculations

The following calculations now have the greatest potential payoff.

1. **Test the \(3\)-adic slopes of the new Padé realizations.**

   In particular:
   \[
   (1-4g^2)^{2n},
   \]
   the quadratic candidate
   \[
   (1-10g^2+23g^4)^n,
   \]
   and the boundary-cancelled Jacobi/Racah families.

   A positive linear \(v_3\)-slope would be qualitatively new.

2. **Prove the new cyclotomic Apéry recurrence.**

   The recurrence verified on 61 terms should be derived by creative telescoping or a Picard–Fuchs argument.

3. **Prove sharp companion denominator theorems.**

   Empirical \(L_{2n}^2\)- and \(L_{12n}^2\)-type bounds should be replaced by exact local \(p\)-adic valuations.

4. **Construct a hybrid radial–Jacobi numerator.**

   The appropriate finite-dimensional search should optimize
   \[
   \log(\text{real weighted norm})
   +
   \log(\text{primitive arithmetic multiplier})
   \]
   directly.

5. **Compare EMN and modular Catalan \(p\)-adic realizations.**

   If both are realizations of the same extension class, an explicit comparison of their syntomic regulators may expose an arithmetic normalizer unavailable in either realization separately.

---

# 26. Principal new results from this investigation

The most useful outcomes can be summarized succinctly.

### Exact or rigorously derived

\[
\int_0^1
{}_2F_1\!\left(\frac14,\frac34;1;t\right)^2dt
=
\frac{16G}{\pi^2}.
\]

The associated symmetric-square moment recurrence.

The moving-EMN moment recurrence and Picard–Fuchs equation

\[
z(z-2)H''+(z-1)H'=\frac1{2(z-1)}.
\]

The higher-pole recurrence

\[
I_{m,t}
=
-\frac1{
2^{t+1}t(t-1)
\binom{2(m-t)}{m-t}}
+
\frac{(m-1)^2}{4t(t-1)}I_{m-2,t-2}.
\]

The Catalan coefficient

\[
b_{m,t}
=
4^{-m}
\binom m{m/2}\binom mt.
\]

The central-pole asymptotics

\[
b_{4n,2n}\sim\frac1{2\pi n},
\qquad
I_{4n,2n}=64^{-n+o(n)}.
\]

The radial weighted problem

\[
W(g)=
\frac{\min(g,1-g)^4}{16g^2}.
\]

The free-integration estimate

\[
256^nL_n
\ll
\left(\frac{16}{27}\right)^n
\]

for the \((1-4g^2)^{2n}\) radial family.

The Jacobi finite-difference cancellation mechanism.

The Racah quadratic-lattice cancellation weights.

The explicit pure-\(G\) simple-pole polynomial

\[
1+569x^2y^2-2800x^4y^4.
\]

### Computationally discovered and strongly verified

The new second-order cyclotomic Catalan Apéry recurrence.

Its common primary/companion recurrence through 61 exact terms.

The denominator laws such as

\[
\operatorname{den}(V_n)\mid d_{2n}^2
\]

over the tested range.

The strong \(2\)-adic and apparently zero \(3\)-adic slopes of the natural EMN approximants.

The degree-two radial candidate

\[
1-10g^2+23g^4.
\]

---

# 27. Conclusion

The search for a qualitatively different Catalan realization succeeded.

The EMN mixed motive provides a genuinely different two-dimensional geometric realization of the Catalan extension, and the direct cyclotomic experiments provide an independent Apéry-type realization.

Neither yet proves

\[
G\notin\mathbf Q.
\]

However, the work identifies a significantly more precise obstruction than before.

For the new EMN Padé realizations it is possible to achieve

\[
\boxed{
\text{integral \(G\)-coefficient}
+
\text{geometrically decaying real form}.
}
\]

What remains is to make the **rational companion** equally efficient.

Thus the Catalan problem has been reduced from

\[
\text{"find a better approximation to \(G\)"}
\]

to

\[
\boxed{
\text{"find a low-height integral lattice for the EMN extension."}
}
\]

The Jacobi/Racah boundary-cancellation identities and the radial weighted-Chebyshev construction give two complementary mechanisms for attacking precisely that lattice problem.

The most promising possible breakthrough would now be a numerator for which those mechanisms interact with a nonzero \(3\)-adic gain. Such a construction would supply a genuinely new adelic Catalan realization and could evade the obstruction that defeated all of the earlier modular hosts.

---

## Reference

P. Eskandari, V. Kumar Murty and Y. Nemoto, *Mixed motives and linear forms in the Catalan constant*, arXiv:2510.20648 (2025).

Internal project references especially relevant to this report include `CATALAN_OBSTRUCTION.md`, `CATALAN_MU4.md`, `CATALAN_EXPLICIT.md`, `MULTI_PRIME_LATTICE.md`, `ADELIC_HOLONOMY.md`, `ONE_CLASS_TWO_WORLDS.md`, and `POSITIVITY_PROGRAM.md`.