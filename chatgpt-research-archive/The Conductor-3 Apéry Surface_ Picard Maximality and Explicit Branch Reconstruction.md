# The Conductor-3 Apéry Surface: Picard Maximality and Explicit Branch Reconstruction

## 1. Main results

Let \(W\) denote the rank-four symplectic local system obtained as the spin/Kuga–Satake lift of the conductor-\(3\) rank-five Apéry hypergeometric system

\[
V=
{}_5F_4\!\left(
\begin{matrix}
\frac12,\frac16,\frac16,\frac56,\frac56\\
1,1,\frac23,\frac43
\end{matrix};z
\right),
\]

so that

\[
V\simeq\Lambda^2_0W.
\]

After the sixfold adapted cover

\[
E:\qquad v^2=u^3-1,\qquad z=u^3,
\]

the finite parts of the local monodromies disappear and \(W\) is realized by the first cohomology of a semistable genus-two fibration

\[
f:X\longrightarrow E.
\]

The calculations below prove:

> **Theorem A — Geometry of the Apéry surface.**
> The associated smooth surface \(X\) has
> \[
> \boxed{p_g(X)=q(X)=1,\qquad K_X^2=2,\qquad c_2(X)=10.}
> \]
> Its genus-two Albanese fibration has only nonseparating degeneration:
> \[
> \boxed{\delta_0=10,\qquad\delta_1=0.}
> \]

The principal new result is:

> **Theorem B — Picard maximality.**
> \[
> \boxed{\rho(X)=10=h^{1,1}(X).}
> \]
> Thus the conductor-\(3\) Apéry surface is Picard-maximal.

This corrects the earlier conjectural value \(\rho=7\).

Finally:

> **Theorem C — Corrected transcendental decomposition.**
> Put
> \[
> H_{\mathrm{par}}=
> H^1(E,j_*W).
> \]
> Then
> \[
> \boxed{\dim H_{\mathrm{par}}=5}
> \]
> but
> \[
> \boxed{
> H_{\mathrm{par}}
> \simeq
> T(X)\oplus\mathbf Q(-1)^3
> }
> \]
> as rational polarized Hodge structures. In particular
> \[
> \boxed{\operatorname{rank}T(X)=2.}
> \]

Thus the earlier hoped-for equality
\[
T(X)\stackrel?=H^1(E,j_*W)
\]
is false: the parabolic cohomology contains three additional algebraic classes.

This is stronger geometrically than the original conjecture.

---

# 2. The adapted cover

The local spectra of \(W\) on

\[
\mathbf P^1-\{0,1,\infty\}
\]

are

\[
T_0=(\omega J_2,\omega^{-1}J_2),
\]

\[
T_1=(-1,-1,1,1),
\]

\[
T_\infty=(-J_2,\zeta_6,\zeta_6^{-1}),
\]

where \(\omega=e^{2\pi i/3}\).

Consequently the minimal common cover eliminating the finite semisimple parts has ramification indices

\[
3,\qquad2,\qquad6.
\]

It is the elliptic curve

\[
\boxed{
E:v^2=u^3-1,\qquad z=u^3.
}
\tag{2.1}
\]

The deck group is cyclic of order six. A convenient generator is

\[
\sigma:(u,v)\longmapsto(\zeta_3u,-v).
\]

The points over the three singularities are:

- two points of ramification \(3\) over \(z=0\);
- three points of ramification \(2\) over \(z=1\);
- one point of ramification \(6\) over \(z=\infty\).

---

# 3. Integral symplectic monodromy

There is an integral lattice on which the monodromies become

\[
M_{0,\mathbf Z}=
\begin{pmatrix}
-1&2&-1&-1\\
1&-4&2&2\\
1&8&-4&-3\\
-6&-19&9&7
\end{pmatrix},
\]

\[
M_{\infty,\mathbf Z}=
\begin{pmatrix}
-8&-5&2&1\\
14&9&-4&-2\\
-22&-15&5&3\\
45&32&-12&-7
\end{pmatrix},
\]

and

\[
M_{1,\mathbf Z}=
\begin{pmatrix}
2&1&0&0\\
-3&-2&0&0\\
0&0&1&0\\
0&0&-3&-1
\end{pmatrix}.
\]

They preserve the unimodular alternating form

\[
J_{\mathbf Z}=
\begin{pmatrix}
0&0&0&1\\
0&0&1&1\\
0&-1&0&0\\
-1&-1&0&0
\end{pmatrix},
\qquad
\det J_{\mathbf Z}=1.
\tag{3.1}
\]

Thus the weight-one VHS has a genuine principal integral polarization.

---

# 4. Stable reduction at \(z=0\)

After the cubic ramification,

\[
N_0=M_{0,\mathbf Z}^{\,3}-I
\]

is square-zero of rank \(2\).

The monodromy pairing

\[
S_0=J_{\mathbf Z}N_0
\]

is positive semidefinite of rank \(2\), with Smith form

\[
\operatorname{SNF}(S_0)
=
\operatorname{diag}(1,3,0,0).
\]

On the primitive rank-two quotient its quadratic form is equivalent over
\(\mathbf Z\) to

\[
\boxed{
A_2=
\begin{pmatrix}
2&-1\\
-1&2
\end{pmatrix}.
}
\tag{4.1}
\]

This is precisely the cycle pairing of the theta graph: two rational components joined at three nodes.

Hence each of the two points above \(z=0\) contributes

\[
\delta_0=3
\]

and has two irreducible components.

Therefore

\[
\boxed{\delta_0(z=0)=6.}
\tag{4.2}
\]

The two fibers also contribute

\[
(2-1)+(2-1)=2
\]

extra vertical divisor classes.

---

# 5. Stable reduction at infinity

At infinity,

\[
N_\infty=M_{\infty,\mathbf Z}^{\,6}-I
=
\begin{pmatrix}
0&0&0&0\\
0&0&0&0\\
8&4&0&0\\
-16&-8&0&0
\end{pmatrix}.
\]

It has rank one and Smith normal form

\[
\boxed{\operatorname{diag}(4,0,0,0).}
\]

Thus it is four times a primitive Picard–Lefschetz operator.

Geometrically, the stable node has thickness \(4\). Its contribution to the boundary is therefore

\[
\boxed{\delta_0(\infty)=4.}
\tag{5.1}
\]

Resolving the local model corresponding to thickness \(4\) inserts three vertical rational components, so this fiber contributes three independent additional component classes.

At \(z=1\), the quadratic ramification kills \(M_1^2=I\), and there is no unipotent degeneration.

Hence altogether

\[
\boxed{\delta_0=6+4=10.}
\tag{5.2}
\]

---

# 6. Hodge degree

For the rank-five system \(V\), take hypergeometric parameters

\[
\alpha=
\left(0,0,0,\frac13,\frac23\right),
\]

\[
\beta=
\left(
\frac16,\frac16,\frac12,\frac56,\frac56
\right).
\]

Applying Fougeron's explicit formula for the parabolic degrees of hypergeometric Hodge bundles gives

\[
\boxed{
\deg_{\mathrm{par}}F^{2,0}_V=\frac16.
}
\tag{6.1}
\]

The relevant parabolic-degree formulas are established for hypergeometric VHSs in Fougeron's work.

The adapted cover has degree \(6\) and removes all parabolic weights, hence

\[
\deg F^{2,0}_{V|E}=1.
\tag{6.2}
\]

Since

\[
V=\Lambda^2_0W
\]

and \(W\) has Hodge numbers \((2,2)\),

\[
F^{2,0}_V
=
\det F^{1,0}_W.
\]

Therefore the Hodge bundle of the genus-two fibration satisfies

\[
\boxed{
\lambda
:=
\deg\det f_*\omega_{X/E}
=1.
}
\tag{6.3}
\]

---

# 7. There are no separating nodes

For a semistable genus-two family, the relation in
\(\operatorname{Pic}(\overline{\mathcal M}_2)\) is

\[
\boxed{
10\lambda=\delta_0+2\delta_1.
}
\tag{7.1}
\]

This is Mumford's genus-two relation.

We have independently computed

\[
\lambda=1,\qquad\delta_0=10.
\]

Hence

\[
10=10+2\delta_1,
\]

so

\[
\boxed{\delta_1=0.}
\tag{7.2}
\]

Thus every degeneration is nonseparating.

This is a useful consistency check: the Hodge calculation and the integral Picard–Lefschetz calculation close exactly.

---

# 8. Numerical invariants of the total surface

Since the base \(E\) has genus one,

\[
\chi(\mathcal O_X)
=
\deg f_*\omega_{X/E}
=1.
\]

The monodromy of \(W\) is irreducible and Zariski-dense in \(\operatorname{Sp}_4\), so \(W\) has no invariant part. Therefore

\[
q(X)=g(E)=1.
\]

Thus

\[
p_g
=
\chi+q-1
=
1.
\]

Moreover, for a semistable genus-two fibration over an elliptic base,

\[
e(X)=\delta=10.
\]

Hence

\[
c_2(X)=10.
\]

Noether's formula gives

\[
K_X^2
=
12\chi(\mathcal O_X)-c_2(X)
=
12-10=2.
\]

Therefore

\[
\boxed{
(p_g,q,K^2,c_2)
=
(1,1,2,10).
}
\tag{8.1}
\]

This places \(X\) exactly in the Bombieri–Catanese–Horikawa class of
\(p_g=q=1,\ K^2=2\) surfaces. Lewis–Lyons describe this moduli space through double covers of the symmetric square of the Albanese elliptic curve.

The Betti numbers follow:

\[
b_1=2,
\]

and from

\[
e(X)=2-2b_1+b_2
\]

we obtain

\[
\boxed{b_2=12.}
\]

Since \(p_g=1\),

\[
\boxed{h^{1,1}=10.}
\tag{8.2}
\]

Thus \(10\) is the absolute Hodge-theoretic upper bound for the Picard number.

---

# 9. The first seven algebraic classes

For a semistable curve fibration, the Leray filtration together with the vertical components gives, over \(\mathbf Q\),

\[
H^2(X)
\simeq
\mathbf Q(-1)^{\,2+\sum_s(m_s-1)}
\oplus
H^1(E,j_*W),
\tag{9.1}
\]

after splitting in the semisimple category of polarized rational Hodge structures.

There are two universal algebraic classes, represented for example by a fiber and a relative polarization.

The component excesses are

\[
1+1+3=5.
\]

Hence

\[
\boxed{
H^2(X,\mathbf Q)
\simeq
\mathbf Q(-1)^7
\oplus H_{\mathrm{par}},
}
\tag{9.2}
\]

where

\[
H_{\mathrm{par}}:=H^1(E,j_*W).
\]

Since \(b_2=12\),

\[
\boxed{\dim H_{\mathrm{par}}=5.}
\tag{9.3}
\]

This explains where the earlier number \(7\) came from: it is the dimension of the immediately visible algebraic subspace.

It is **not** the Picard number.

---

# 10. The crucial \(C_6\)-decomposition

The cover \(E\to\mathbf P^1\) is cyclic of degree six.

Let

\[
\chi_k(\sigma)=\zeta_6^k,
\qquad0\le k<6.
\]

By the projection formula,

\[
H_{\mathrm{par},\mathbf C}
=
\bigoplus_{k=0}^{5}H_k,
\]

with

\[
H_k=
H^1_{\mathrm{par}}
\left(
\mathbf P^1,
W\otimes L_{\chi_k}
\right).
\tag{10.1}
\]

For a rank-\(n\) local system on the thrice-punctured sphere, the parabolic Euler-characteristic formula expresses the dimension through the local fixed spaces.

Applying it to the three twisted local monodromies gives

\[
\boxed{
(\dim H_0,\ldots,\dim H_5)
=
(2,0,1,1,1,0).
}
\tag{10.2}
\]

Thus, over \(\mathbf Q\), the \(C_6\)-representation decomposes as

\[
\boxed{
H_{\mathrm{par}}
=
H_{\mathrm{triv}}
\oplus
H_{\mathrm{sgn}}
\oplus
H_{\Phi_3},
}
\tag{10.3}
\]

with rational dimensions

\[
\boxed{2+1+2=5.}
\]

Here:

- \(H_{\mathrm{triv}}\) has rank \(2\);
- \(H_{\mathrm{sgn}}\) has rank \(1\);
- \(H_{\Phi_3}\) is the rank-\(2\) rational piece combining \(k=2,4\).

---

# 11. Proof of Picard maximality

From (9.2), the Hodge numbers of the parabolic part must be

\[
\boxed{
h(H_{\mathrm{par}})=(1,3,1).
}
\tag{11.1}
\]

Now use the rational \(C_6\)-decomposition.

A rank-one rational polarized Hodge structure of weight two can only be of type

\[
(1,1).
\]

Therefore

\[
H_{\mathrm{sgn}}
\simeq
\mathbf Q(-1).
\]

Of the two remaining rank-two rational pieces, exactly one contains the unique \(H^{2,0}\)-line.

Because it is defined over \(\mathbf Q\), complex conjugation then forces it also to contain the \(H^{0,2}\)-line. Since its rank is exactly two, its entire Hodge decomposition is

\[
(1,0,1).
\]

Call this piece \(T\).

The other rank-two summand is therefore entirely of type

\[
(0,2,0),
\]

i.e. two rational \((1,1)\)-classes.

By the Lefschetz \((1,1)\)-theorem, all three rational \((1,1)\)-directions are algebraic.

The seven already-visible algebraic classes therefore acquire three more:

\[
\rho(X)\ge7+3=10.
\]

But

\[
\rho(X)\le h^{1,1}(X)=10.
\]

Hence

\[
\boxed{\rho(X)=10.}
\tag{11.2}
\]

This proves Theorem B.

The surface is Picard-maximal.

---

# 12. The corrected transcendental statement

Since

\[
b_2=12,\qquad\rho=10,
\]

the transcendental lattice has

\[
\boxed{\operatorname{rank}T(X)=2.}
\tag{12.1}
\]

The preceding argument actually proves the rational Hodge decomposition

\[
\boxed{
H^1(E,j_*W)
\simeq
T(X)\oplus\mathbf Q(-1)^3.
}
\tag{12.2}
\]

Thus the rank-five parabolic cohomology is itself an Apéry-type object with a particularly clean split:

\[
\boxed{
5=2_{\rm transcendental}+3_{\rm algebraic}.
}
\]

This is the corrected form of the earlier proposed identification.

A rank-two weight-two transcendental Hodge structure is the smallest possible nontrivial one for a surface with \(p_g=1\). Consequently the associated K3-type Hodge structure is of singular-K3 type: any K3 surface carrying this same transcendental structure has Picard number \(20\).

This is considerably more rigid than the initially anticipated rank-five transcendental geometry.

---

# 13. Interpretation for the Apéry theory

The original conductor-\(3\) hypergeometric object had rank five and Hodge numbers

\[
(1,3,1).
\]

The present calculation explains the five dimensions geometrically:

\[
\boxed{
V_{\rm rank\,5}
\rightsquigarrow
T(X)_{\rm rank\,2}
\oplus
\mathbf Q(-1)^3.
}
\]

Thus three directions in the apparently transcendental rank-five Apéry system become divisor classes after the Kuga–Satake/genus-two realization and the sixfold cyclotomic base change.

Only a rank-two core remains genuinely transcendental.

This suggests a refinement of the general Apéry philosophy:

> **Arithmetic rank can be strictly smaller than differential-equation rank.**

The relevant quantity for irrationality and period independence may therefore be the rank of the non-Tate quotient after all geometric realizations and finite covers, rather than the raw rank of the Picard–Fuchs operator.

For conductor \(3\), that effective transcendental rank is only

\[
\boxed{2.}
\]

---

# 14. Explicit branch reconstruction

The Bombieri–Catanese–Horikawa description writes a
\(p_g=q=1,\ K^2=2\) surface as a double cover

\[
X\longrightarrow E^{(2)}
\]

of the symmetric square of its Albanese elliptic curve, with branch divisor in a fixed six-dimensional projective linear system. Lewis–Lyons give explicit polynomial realizations of this construction.

Our monodromy computation independently shows that the Albanese curve is equianharmonic:

\[
\boxed{j(E)=0.}
\]

The mod-\(2\) representation

\[
\operatorname{Sp}_4(\mathbf F_2)\simeq S_6
\]

acts on the six Weierstrass labels with local cycle types

\[
(6),\qquad(2,2,2),\qquad(3,1,1,1).
\]

A compatible Belyi cover is

\[
\boxed{
\phi(x)=-\frac{x^6}{2x^3+1}.
}
\tag{14.1}
\]

After pullback to \(E\), a normalization of this label cover can be written

\[
\boxed{
C:\quad y^6=-(2x^3+1),
}
\tag{14.2}
\]

with

\[
\boxed{
u=\frac{x^2}{y^2},
\qquad
v=\frac{x^3+1}{y^3}.
}
\tag{14.3}
\]

This correctly reconstructs the permutation cover of the six Weierstrass labels.

It does **not**, by itself, reconstruct their cross-ratios. Thus (14.2) is not yet the equation of the genus-two family.

---

# 15. A striking explicit candidate from Lyons's high-Picard family

Christopher Lyons has recently constructed a one-parameter family of
\(p_g=q=1,\ K^2=2\) surfaces with Picard number at least \(8\), using explicit branch equations.

His auxiliary elliptic curve is

\[
E_\tau:
y^2=
x^3+
(1-6\tau^2-3\tau^4)x^2
+
16\tau^6x.
\tag{15.1}
\]

The excluded parameters are

\[
\tau\in\{0,\pm1,\pm1/3\},
\]

so the specialization used below is legitimate.

Write

\[
A=1-6\tau^2-3\tau^4,
\qquad
B=16\tau^6.
\]

The quotient by the rational two-torsion point \((0,0)\) has model

\[
\widehat E_\tau:
y^2=x^3-2Ax^2+(A^2-4B)x.
\]

Its \(c_4\)-invariant is proportional to

\[
A^2+12B.
\]

A direct factorization gives

\[
\boxed{
A^2+192\tau^6
=
(3\tau^2+1)
(3\tau^6+75\tau^4-15\tau^2+1).
}
\tag{15.2}
\]

Consequently

\[
\boxed{\tau^2=-\frac13}
\tag{15.3}
\]

makes the Albanese quotient equianharmonic:

\[
j(\widehat E_\tau)=0.
\]

This exactly matches the Albanese curve forced by the Apéry monodromy.

---

# 16. Candidate branch coefficients

Lyons's branch section is written

\[
\mu(\tau)
=
c_0\lambda_0+
c_1\lambda_1+
c_2\lambda_2+
c_3\lambda_3,
\]

with explicit polynomial coefficients in \(\tau\).

Substitution of

\[
\tau^2=-\frac13
\]

and removal of a common scalar gives

\[
\boxed{
\begin{aligned}
c_0&=24(1-45\tau),\\
c_1&=72(7+9\tau),\\
c_2&=8(35-39\tau),\\
c_3&=27(9\tau-13).
\end{aligned}}
\tag{16.1}
\]

Thus we obtain the completely explicit candidate

\[
\boxed{
\begin{aligned}
B_{\rm cand}=Z\big(&
24(1-45\tau)\lambda_0
+
72(7+9\tau)\lambda_1\\
&+
8(35-39\tau)\lambda_2
+
27(9\tau-13)\lambda_3
\big),
\\
&\hspace{25mm}\tau^2=-\frac13.
\end{aligned}}
\tag{16.2}
\]

This candidate lies in a published explicit family specifically constructed to have unusually large Picard number.

It also has exactly the equianharmonic Albanese curve required by our independent Apéry calculation.

These coincidences make (16.2) a serious candidate for the missing explicit branch divisor.

---

# 17. What is and is not proved about the branch equation

It is important not to conflate the preceding evidence with an identification theorem.

We have proved:

1. the Apéry surface has \(j(\operatorname{Alb}X)=0\);
2. it is \(p_g=q=1,\ K^2=2\);
3. it is Picard-maximal, \(\rho=10\);
4. its genus-two monodromy and complete stable degeneration data are known;
5. its six-Weierstrass permutation cover is known;
6. Lyons's family contains an explicit nonexcluded \(j=0\) specialization;
7. its branch coefficients are explicitly (16.1).

We have **not yet proved**

\[
\boxed{
X_{\rm Apéry}\cong X_{\tau},
\qquad
\tau^2=-1/3.
}
\tag{17.1}
\]

Equality of the Albanese elliptic curve and a large Picard number is insufficient to identify two surfaces.

One further exact comparison is required. Any one of the following would finish it:

\[
\boxed{
\begin{array}{ll}
\text{(i)}&
T(X_\tau)\cong T(X_{\rm Apéry});\\[1mm]
\text{(ii)}&
R^1f_{\tau,*}\mathbf Z
\cong W_{\rm Apéry};\\[1mm]
\text{(iii)}&
\text{the genus-two Picard--Fuchs operator of }X_\tau
\text{ equals }\mathscr T;\\[1mm]
\text{(iv)}&
\text{their period points/Igusa invariants agree.}
\end{array}}
\tag{17.2}
\]

Until such a comparison is made, (16.2) should be labelled **candidate**, not theorem.

---

# 18. Final theorem ledger

### Proved

\[
\boxed{\deg F^{2,0}_{V|E}=1}
\]

\[
\boxed{\lambda=1}
\]

\[
\boxed{\delta_0=10,\quad\delta_1=0}
\]

\[
\boxed{p_g=q=1,\quad K^2=2,\quad c_2=10}
\]

\[
\boxed{b_2=12,\quad h^{1,1}=10}
\]

\[
\boxed{\dim H^1(E,j_*W)=5}
\]

\[
\boxed{
H^1(E,j_*W)
\simeq
T(X)\oplus\mathbf Q(-1)^3
}
\]

\[
\boxed{\operatorname{rank}T(X)=2}
\]

and, most importantly,

\[
\boxed{\rho(X)=10.}
\]

### Disproved

The earlier conjectures

\[
\rho(X)=7
\]

and

\[
T(X)=H^1(E,j_*W)
\]

are false.

Their corrected versions are the Picard-maximality and decomposition theorems above.

### Explicit but not yet identified

The Lyons specialization

\[
\boxed{\tau^2=-1/3}
\]

gives an equianharmonic \(p_g=q=1,\ K^2=2\) surface with explicit candidate branch divisor (16.2).

The final equality with the Apéry surface remains to be established by an exact transcendental-lattice or Picard–Fuchs comparison.

---

# 19. Conceptual conclusion

The conductor-\(3\) Apéry construction now has the structure

\[
\boxed{
\begin{array}{ccccc}
L(2,\chi_{-3})
&\longleftarrow&
\text{Apéry extension}
&\longleftarrow&
V_5
\\
&&&&\rotatebox{90}{$\simeq$}\\[-1mm]
&&&&\Lambda^2_0W_4
\\[1mm]
&&&&\downarrow\\
&&&&
\text{genus-two VHS over }\mathbf P^1
\\
&&&&\downarrow\text{ sixfold cover}\\
&&&&
f:X\to E_{j=0}
\\
&&&&\downarrow\\
&&&&
(p_g,q,K^2,\rho)=(1,1,2,10).
\end{array}}
\]

After passing to the surface, the apparent rank-five K3-type object splits as

\[
\boxed{
V_{\rm geometric}
=
T(X)_{\operatorname{rank}2}
\oplus
\mathbf Q(-1)^3.
}
\]

Thus the deepest arithmetic content of the conductor-\(3\) Apéry system is carried by a rank-two transcendental Hodge structure on a Picard-maximal surface.

This suggests a revised principle for Apéry-sequence theory:

\[
\boxed{
\textbf{The effective arithmetic dimension of an Apéry realization is the rank of its non-Tate quotient, not necessarily the order of its differential equation.}
}
\]

For the conductor-\(3\) cellular construction, that dimension is exactly

\[
\boxed{2.}
\]

That is the strongest geometric structural result obtained from this branch of the theory so far.