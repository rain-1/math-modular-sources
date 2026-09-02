# Elliptic hosts: fold periods, denominator types, and closed forms

Working directory: `/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/elliptic`

Scripts / logs in that directory:
`task1.gp`,`task1.out` (Task 1) · `periods.py`,`master.py`,`consts.py`,`check.py`,`D.pkl`,`fp.pkl`,`fp.out` (Tasks 2, 7) ·
`ell.gp`,`ell.out`,`lval.gp`,`lval.out`,`gpt.gp` (Task 3) ·
`pslq1.py`,`pslq1.out`,`scan.py`,`scan.out`,`scan2.out`,`gscan.out`,`gscan2.out`,`fpscan.out`,`final_scan.out` (Task 4) ·
`boyd.gp`,`boyd.out`,`mahler.out` (Task 5) · `tables.out` (all numerics at 60 digits).

Throughout, for a host
$$H(x)=\bigl((1-a_0x)(1-ax)(1-bx)\bigr)^{-1/2},\qquad P(x)=(1-a_0x)(1-ax)(1-bx),$$
write $r_1<r_2<r_3$ for the roots of $P$ (the branch points), and
$$\Omega_{ij}=\int_{r_i}^{r_j}\frac{dx}{\sqrt{|P|}},\qquad \eta_{ij}=\int_{r_i}^{r_j}\frac{x\,dx}{\sqrt{|P|}} .$$
For the $E_m,F_m$ hosts $a_0=1$, so $r_3=1$: the puncture **is** a branch point.
For $G$, $(a_0,a,b)=(4,8,12)$ and $r_3=1/4<1$: the puncture $x=1$ is a regular point.

---

## 0. Branch structure of the continued host (needed for everything below)

Continuing $H$ from $x=0$ along the real axis in the upper half plane
($x=u+i0^+$), each linear factor $1-\alpha u$ that has turned negative
contributes $(\,\cdot\,)^{-1/2}=i\,|\cdot|^{-1/2}$. Hence

$$H(u+i0^+) \;=\; i^{\,\#\{\alpha:\ \alpha u>1\}}\ |P(u)|^{-1/2}.$$

**[verified numerically, 20 digits, `check.py`]** — checked at $u=0.15,0.2$ for $E_1$
and $u=0.1,0.2$ for $G$ against direct evaluation of $1/\sqrt{P(u+10^{-25}i)}$.

Consequences for the fold periods $c_i[k]=\int_0^{\delta_i}Hk\,dt$:

| host type | $c$ at $r_1$ | $c$ at $r_2$ | $c$ at $r_3$ |
|---|---|---|---|
| $E_m,F_m$ | $A[k]$ (real) | $A[k]+i\,B[k]$ | (puncture) |
| $G$ | $A[k]$ (real) | $A[k]+i\,B[k]$ | $A[k]-C[k]+i\,B[k]$ |

with the **real** integrals
$$A[k]=\int_0^{r_1}\frac{k\,dt}{\sqrt{|P|}},\quad B[k]=\int_{r_1}^{r_2}\frac{k\,dt}{\sqrt{|P|}},\quad C[k]=\int_{r_2}^{r_3}\frac{k\,dt}{\sqrt{|P|}}.$$

So for $E_m,F_m$ **the difference period is purely imaginary**:
$$\Delta[k]\;=\;c_1[k]-c_2[k]\;=\;\pm\, i\,B[k],$$
the sign being $+$ when $\delta_1=1/9$ (resp. $1/25$) is the *outer* fold
($m\ge3$ for $E$; never for $F$) and $-$ otherwise.
Explicitly: $\Delta=-iB$ for $E_1,E_2,F_1,F_2,F_3$ and $\Delta=+iB$ for $E_3,E_4,E_5,E_6,E_{12}$.

Numerically all integrals were computed with **analytically regularised, cancellation-free
integrands** (endpoint substitutions $t=r-s^2$ and $t=r_i+(r_j-r_i)\frac{1-\cos\varphi}{2}$
which cancel the inverse-square-root exactly), at `mp.dps=90`–`110`.
Agreement between `dps=90` and `dps=140` is **>75 digits**
(`A_B(E_1)`, `B_D(E_1)` checked); 60 digits are quoted below.

---

## 1. Integrality and denominator types  [verified exact, $n\le 300$, PARI]

Output: `task1.out`.

For **all 11 hosts** ($E_m$, $m\in\{1,2,3,4,5,6,12\}$; $F_m$, $m\in\{1,2,3\}$; $G$):

* $H\in\mathbb Z[[x]]$: **YES**, all coefficients $n=0,\dots,300$ are integers.
* $H_B=H\!\int_0^x H/(1-t)\,dt$: type $[1..n]$, **max excess $=1$** (i.e. $\operatorname{lcm}(1..n)\,[x^n]H_B\in\mathbb Z$ for all $n\le300$).
* $H_D=H\!\int_0^x H\log(1-t)/(1-t)\,dt$: type $[1..n][1..n/2]$, **max excess $=1$**.
* $H_L=H\!\int_0^x H\log(1-t)\,dt$: type $[1..n][1..n/2]$, **max excess $=1$**.

The second factor is genuinely needed: the max excess of $H_D$ (resp. $H_L$) over
$\operatorname{lcm}(1..n)$ *alone* is astronomically large, e.g. for $E_1$ at $n=300$

* $H_D$: excess $= 8356221165259964192329704125048709983866469634630441370296000$,
* $H_L$: excess $= 253556694314680143777807792007094321572562103286032000$,

and comparable numbers ($10^{50}$–$10^{62}$, attained at $n=287$–$300$) for every other host.

First coefficients of $H$ (a check that the hosts are the intended ones):

| host | $[x^0..x^7]H$ |
|---|---|
| $E_1$ | 1, 7, 49, 361, 2779, 22099, 179875, 1488823 |
| $E_2$ | 1, 9, 77, 657, 5613, 48033, 411725, 3534993 |
| $E_3$ | 1, 11, 117, 1253, 13543, 147687, 1623975, 17994747 |
| $E_4$ | 1, 13, 169, 2269, 31369, 444301, 6416425, 94099837 |
| $E_5$ | 1, 15, 233, 3825, 65571, 1160235, 21002555, 386496975 |
| $E_6$ | 1, 17, 309, 6041, 124309, 2647209, 57689205, 1277437593 |
| $E_{12}$ | 1, 29, 1017, 39917, 1660633, 71377053, 3131105145, 139278861453 |
| $F_1$ | 1, 15, 273, 5585, 121131, 2713275, 62012355, 1436993775 |
| $F_2$ | 1, 17, 317, 6441, 138349, 3077257, 70011229, 1617522889 |
| $F_3$ | 1, 19, 373, 7693, 165079, 3649999, 82512871, 1896108787 |
| $G$ | 1, 12, 128, 1344, 14208, 152064, 1648640, 18087936 |

---

## 2. Fold periods  [numerical, 60 digits quoted, $\ge75$ digits accurate]

Full 60-digit table: `tables.out`. Selected values ($\Omega=\Omega_{12}$, $\eta=\eta_{12}$):

### $E_1$: $P=(1-x)(1-9x)(1-4x)$, folds $1/9,1/4$
```
A_B = 0.303167065196972484877931850241413804016090345628367486705413
B_B = 0.710647653194800437952723541739195475339440139413410499327335
A_D = -0.0250018730365440417398577923074866771386140226588845530526433
B_D = -0.146633077532013280740455612601978496558655464756834105539644
A_L = -0.0227131337940616272218873665409640362200078954271850879048227
B_L = -0.117421895661402275231267555618887979780352728954094897162198
Omega_12 = 0.579196055657041066671237163795353181379731945102973260906671
eta_12   = 0.105430953527174108036088135969222864486771852160699594688448
```

### $G$: $P=(1-4x)(1-8x)(1-12x)$, folds $1/12,1/8,1/4$
```
A_B = 0.287243980831620684807362925163284254123261959219818001389726
B_B = 0.470962632642101546025683046304949890851006989489441175170776
C_B = 0.656927507758595418866018732182621265823645172016784906741637
A_D = -0.0183293368302849223996314767723999923419067857604626204722676
B_D = -0.0523909042022578426679653477032971431849292605462964003044072
C_D = -0.130788650567804406886249785965979956047678772701627632497449
A_L = -0.0170460611115184257613556104215434949995598347642308453428555
B_L = -0.0467677262025409986300893768766522500453850493313780610634925
C_L = -0.105753374049536668583228271700397621663280993805804078394325
Omega_01 = 0.26956445593745540442983437485004025360801377706160317004562
Omega_12 = 0.421437588703149010717800914449769247375200223535272261029987
Omega_23 = 0.53912891187491080885966874970008050721602755412320634009124
```

Summary of the fold-difference magnitudes $|\Delta[k]|=B[k]$:

| host | $B_B$ | $B_D$ | $B_L$ | $\Omega_{12}$ | $\eta_{12}$ |
|---|---|---|---|---|---|
| $E_1$ | 0.71064765319480043795 | −0.14663307753201328074 | −0.11742189566140227523 | 0.57919605565704106667 | 0.10543095352717410804 |
| $E_2$ | 0.44704065978930196984 | −0.05618774266128099313 | −0.04953994110183176410 | 0.39424669280395328553 | 0.04654840185671842010 |
| $E_3$ | 0.35250321596853464197 | −0.03613706613820352823 | −0.03258029020319341936 | 0.31817557698225286285 | 0.03095073434122463606 |
| $E_4$ | 0.30020109414399719907 | −0.02747293815836082480 | −0.02497771613265977113 | 0.27399630160999125650 | 0.02382872315666025728 |
| $E_5$ | 0.26587364992736026356 | −0.02262337313441069636 | −0.02064874531559598152 | 0.24425353105583100457 | 0.01973800445050455979 |
| $E_6$ | 0.24113770502975265940 | −0.01950285827850771885 | −0.01783662344572732405 | 0.22248133528246314070 | 0.01706773470156272714 |
| $E_{12}$ | 0.16781097634155283768 | −0.01184534549367453879 | −0.01086155422265438033 | 0.15646546170293364487 | 0.01040776007232284171 |
| $F_1$ | 0.40307829938779551755 | −0.06927055730320794551 | −0.05596469462884747953 | 0.34072451610537582782 | 0.05050814054616305518 |
| $F_2$ | 0.25327968119709795002 | −0.02235212439308351097 | −0.02022586597603641266 | 0.23201001491575805617 | 0.01925508271019577815 |
| $F_3$ | 0.19965062366710191928 | −0.01281424563231942716 | −0.01196928885010824555 | 0.18726391565417960382 | 0.01157136682712991485 |
| $G$ | 0.47096263264210154603 | −0.05239090420225784267 | −0.04676772620254099863 | 0.42143758870314901072 | 0.04421513845331112120 |

---

## 3. Elliptic data  [PARI, 100 digits; quoted to ~30]

Model used: $y^2=P(x)$ with $P=Ax^3+Bx^2+Cx+1$, $x=X/A$, $y=Y/A$, giving the
integral model $Y^2=X^3+BX^2+ACX+A^2$. Note $dx/y = dX/Y = 2\cdot dX/(2Y)$, so PARI's
periods (normalised for $dX/(2Y)$) are **exactly** our interval integrals:

> **[verified numerically, $\ge$60 digits]**  $\Omega_{23}=\omega_1(E)$ and $\Omega_{12}=|\omega_2(E)|$
> for the model $[0,B,0,AC,A^2]$, for every host. (E.g. $E_1$: $\omega=[0.83798706481037411236\ldots,\ -0.57919605565704106667\ldots\,i]$.)

| host | model $[0,B,0,AC,A^2]$ | minimal model | $N$ | $j$ | analytic rank | order of the point $x=0$ |
|---|---|---|---|---|---|---|
| $E_1$ | $[0,49,0,504,1296]$ | $[1,0,1,-19,26]$ | 30 | $702595369/72900$ | 0 | **3** |
| $E_2$ | $[0,89,0,1296,5184]$ | $[1,1,1,-84,261]$ | 42 | $65597103937/63504$ | 0 | **4** |
| $E_3$ | $[0,129,0,2376,11664]$ | $[1,-1,0,-198,1120]$ | 198 | $1180932193/4356$ | 1 | $\infty$ |
| $E_4$ | $[0,169,0,3744,20736]$ | $[1,0,0,-361,2585]$ | 210 | $5203798902289/57153600$ | 0 | **3** |
| $E_5$ | $[0,209,0,5400,32400]$ | $[1,1,0,-572,4956]$ | 6270 | $20753798525641/353816100$ | 1 | $\infty$ |
| $E_6$ | $[0,249,0,7344,46656]$ | $[1,-1,1,-833,9281]$ | 2070 | $87587538121/1904400$ | 1 | $\infty$ |
| $E_{12}$ | $[0,489,0,25056,186624]$ | $[1,-1,1,-3416,75291]$ | 10998 | $6045477024313/215032896$ | 1 | $\infty$ |
| $F_1$ | $[0,129,0,3000,10000]$ | $[1,-1,0,-159,665]$ | 630 | $611960049/122500$ | 1 | $\infty$ |
| $F_2$ | $[0,233,0,6800,40000]$ | $[1,1,1,-706,6719]$ | 3570 | $38920307374369/1274490000$ | 1 | $\infty$ |
| $F_3$ | $[0,337,0,11400,90000]$ | $[1,0,1,-1654,25556]$ | 4290 | $499980107400409/4140922500$ | 1 | $\infty$ |
| $G$ | $[0,176,0,9216,147456]$ | $[0,-1,0,-4,4]$ | **24** | $35152/9$ | 0 | **4** |

$L$-values (minimal model; $L(E,0)=0$ always and $L'(E,0)=\varepsilon\,N\,L(E,2)/(4\pi^2)$
by the functional equation — **[verified numerically, 60 digits]** — so $L(E,2)$ carries no
information beyond $L'(E,0)$):

| host | $L'(E,0)$ | $L(E,2)$ |
|---|---|---|
| $E_1$ | 0.616870938789551946232005208740 | 0.811769617650864831525184556866 |
| $E_2$ | 1.15459137735997967479126480950 | 1.08527239423351236670376345766 |
| $E_3$ | −3.28435857383187368804400533762 | 0.654854946162558176394142445033 |
| $E_4$ | 7.65435423604863322713542854563 | 1.43896091915433458802235999167 |
| $E_5$ | −114.542959054427772845465538612 | 0.721208097478391416925393488848 |
| $E_6$ | −62.5721126625089228725597354082 | 1.19335651887798502502389002771 |
| $E_{12}$ | −338.977337531668040266769777357 | 1.21679295230845750637787929926 |
| $F_1$ | −11.9554569731145308446077604138 | 0.749178607992924669356467924919 |
| $F_2$ | −102.103928920703490320541153867 | 1.12910407422329435759059157213 |
| $F_3$ | −92.1162490128473619987501965983 | 0.847693181043395440846719983313 |
| $G$ | **0.511424067053503722283274426476** | 0.841258870502379931364584375666 |

**Torsion side-remark [verified numerically, 70 digits, PSLQ].** Where the rational point
$(x,y)=(0,1)$ is torsion, the *partial* period $\Omega_{01}=\int_0^{r_1}dx/y$ is a rational
multiple of a full period: $\Omega_{23}=3\,\Omega_{01}$ for $E_1$ and $E_4$ (order 3),
$\Omega_{23}=2\,\Omega_{01}$ for $E_2$ and $G$ (order 4). No such relation for the rank-1
hosts, where the point has infinite order.

---

## 4. Closed forms for the difference periods  [numerical PSLQ, residuals $\sim10^{-91}$]

Write $\;\Delta=(a-1)(b-1)$, $\;f=\dfrac{ab}{\Delta}$, $\;c=\dfrac1f=\dfrac{(a-1)(b-1)}{ab}=(1-r_a)(1-r_b)$,
$\;\Omega=\Omega_{12}$, $\;\eta=\eta_{12}$.  (Here $\{a,b\}=\{9,4m\}$ or $\{25,4m\}$, so $r_a,r_b$ are the two folds.)

### 4a. Three exact identities for $E_m$ and $F_m$ (puncture $=$ branch point)

> **(I)** $\displaystyle B_B=\int_{r_1}^{r_2}\frac{dx}{(1-x)\,y}\;=\;f\,\bigl(\Omega-\eta\bigr)$
>
> **(II)** $\displaystyle B_L=\int_{r_1}^{r_2}\frac{\log(1-x)\,dx}{y}\;=\;\frac{\Omega}{2}\,\log c$
>
> **(III)** $\displaystyle B_D=\int_{r_1}^{r_2}\frac{\log(1-x)\,dx}{(1-x)\,y}\;=\;\frac{(a+b)\,\Omega-2ab\,\eta}{\Delta}\;+\;\frac{f}{2}\,(\Omega-\eta)\,\log c$
>
> and for the near-fold ("$A$") periods
>
> **(IV)** $\displaystyle A_B=f\,(\Omega_{01}-\eta_{01})-\frac{2}{\Delta}.$

All four verified for **all ten** $E_m,F_m$ hosts, residual $\lesssim10^{-90}$ at `mp.dps=110`
(`scan2.out`, `pslq1.out`, `fpscan.out`). Explicit constants:

| host | $a,b$ | $\Delta$ | $f=ab/\Delta$ | $c$ | $(a+b)$ | $2ab$ |
|---|---|---|---|---|---|---|
| $E_1$ | 9,4 | 24 | $3/2$ | $2/3$ | 13 | 72 |
| $E_2$ | 9,8 | 56 | $9/7$ | $7/9$ | 17 | 144 |
| $E_3$ | 9,12 | 88 | $27/22$ | $22/27$ | 21 | 216 |
| $E_4$ | 9,16 | 120 | $6/5$ | $5/6$ | 25 | 288 |
| $E_5$ | 9,20 | 152 | $45/38$ | $38/45$ | 29 | 360 |
| $E_6$ | 9,24 | 184 | $27/23$ | $23/27$ | 33 | 432 |
| $E_{12}$ | 9,48 | 376 | $54/47$ | $47/54$ | 57 | 864 |
| $F_1$ | 25,4 | 72 | $25/18$ | $18/25$ | 29 | 200 |
| $F_2$ | 25,8 | 168 | $25/21$ | $21/25$ | 33 | 400 |
| $F_3$ | 25,12 | 264 | $25/22$ | $22/25$ | 37 | 600 |

In the uniform families: $c(E_m)=\dfrac{2(4m-1)}{9m}$, $c(F_m)=\dfrac{6(4m-1)}{25m}$.

**Mechanism (not proved here, but forced by the shape).** $x=1$ is a Weierstrass point;
translation by the $2$-torsion point $(1,0)$ is an involution $\sigma$ of the curve that
fixes the fold cycle setwise, exchanges $r_a\leftrightarrow r_b$, and satisfies
$(1-x)(1-\sigma x)=c$. Averaging $\log(1-x)$ over $\{x,\sigma x\}$ gives (II) at once.
(I) is the reduction of the second-kind differential $dx/((1-x)y)$ — it has a *double pole
with zero residue* at the Weierstrass point — to the basis $\{dx/y,\ x\,dx/y\}$, and (IV) is
the same reduction with the boundary term $-\tfrac{2}{\Delta}\,\frac{y}{1-x}\big|_{x=0}$.

**Negative results [PSLQ, 65–70 digits, `final_scan.out`].** The *near-fold* periods $A_D$
and $A_L$ are **not** in the span of $\{1,\Omega_{01},\eta_{01},\Omega,\eta,\Omega\log c,\eta\log c\}$
with coefficients $\le10^5$ — they are genuinely new transcendentals. Nor are $B_D,B_L$
expressible through $L'(E,0)$ or $\pi$ alone.

### 4b. $G$ (puncture $\ne$ branch point): no linear relation, one exact bilinear one

For $G$ the differential $dx/((1-x)y)$ is of the **third kind**: $P(1)=(1-4)(1-8)(1-12)=-231$,
so the two points over $x=1$ are $(1,\pm i\sqrt{231})$ and the residues are $\mp1/(i\sqrt{231})$.
Accordingly there is **no** $\mathbb{Q}$-linear relation: PSLQ with `maxcoeff` $10^6$ and tolerance
$10^{-70}$ finds **nothing** for any of $A_B,A_L,A_D,B_B,B_L,B_D,C_B,C_L,C_D$ against
$\{1,\pi,\Omega_{01},\eta_{01},\Omega_{12},\eta_{12},\Omega_{23},\eta_{23},L'(E,0),\Omega\log2,\Omega\log3,\Omega\log7,\Omega\log11,\pi u/\sqrt{231}\}$
in several subsets (`gscan.out`, `gscan2.out`).

Instead the Riemann bilinear relation holds exactly:

> **(V)** $\displaystyle B_B\cdot\Omega_{23}\;-\;C_B\cdot\Omega_{12}\;=\;-\,\frac{2\pi\,u}{\sqrt{231}}$
>
> where $u=0.0555011930409905959917529873356\ldots$ is (minus $i$ times) the elliptic
> logarithm of the point $(-384,\,384\,i\sqrt{231})$ on $Y^2=X^3+176X^2+9216X+147456$,
> i.e. of the point over $x=1$, in PARI's $dX/(2Y)$ normalisation. Equivalently
> $u=|{\rm ellpointtoz}|$, purely imaginary, with $u/\omega_2=0.1316949283327557327848604317\ldots$

**[verified numerically, PSLQ, $\ge70$ digits; the ratio $(B_B\Omega_{23}-C_B\Omega_{12})\sqrt{231}/(\pi u)$
equals $-2.0000\ldots$ to 40 printed digits, PSLQ residual $<10^{-70}$.]**

The analogous bilinear combinations with the $L$- and $D$-kernels
($B_L\Omega_{23}-C_L\Omega_{12}=0.01935461361821922328\ldots$,
$B_D\Omega_{23}-C_D\Omega_{12}=0.02687380235032826712\ldots$)
match **nothing** in $\{\pi u/\sqrt{231},\pi\Omega_{23},\pi\Omega_{12},\Omega_{23},\Omega_{12}\}$
at `maxcoeff` $10^5$, $10^{-65}$.

Also, $\exp(2B_L/\Omega_{12})=0.800960391286662539553879457572\ldots$ is **not** rational
(mpmath `identify` with `maxcoeff` $10^5$: nothing), in sharp contrast with (II), where it
is $c\in\mathbb{Q}$ for all ten $E_m,F_m$.

---

## 5. Mahler measures / Boyd's family  [PARI + numerical, 40 digits]

Boyd's curves for $x+1/x+y+1/y+k=0$ were built with `ellfromeqn` and compared by
conductor and by $a_p$ at the first 12 primes (`boyd.out`).

* **$G$ is isogenous to Boyd's $k=2$ and $k=8$ curves.** All three have conductor $24$
  and the identical $a_p$-vector $[0,-1,-2,0,4,-2,2,-4,-8,6,8,6]$ for $p\le37$.
  They are **not isomorphic**: $j(G)=35152/9$, $j(k{=}2)=2048/3$, $j(k{=}8)=28756228/3$.
  (Isogeny class 24a; $G$'s minimal model is $[0,-1,0,-4,4]$.)
* **No $E_m$ or $F_m$ matches any Boyd curve with $1\le k\le24$** (their conductors
  $30,42,198,210,6270,2070,10998,630,3570,4290$ do not occur in Boyd's list
  $15,24,21,15,120,231,24,195,840,1155,48,663,840,3135,15,4641,1848,6555,240,1785,3432,1311,840$).

Mahler measures, via $m=\frac1\pi\int\operatorname{arccosh}\bigl(|2\cos\theta+k|/2\bigr)d\theta$
over $\{|2\cos\theta+k|\ge2\}$ (`mahler.out`, `mp.dps=60`):

$$m(x+\tfrac1x+y+\tfrac1y+2)=0.5114240670535037222832744264759064461855\ldots = L'(E_{24},0)$$
$$m(x+\tfrac1x+y+\tfrac1y+8)=2.045696268214014889133097705903625784742\ldots = 4\,L'(E_{24},0)$$

**[verified numerically, 40 digits]** — these are Boyd's (proved) $k=2,8$ identities and
$L'(E_{24},0)=0.51142406705350372228327442647590644618552524058474\ldots$ is exactly $G$'s $L'$.

However, PSLQ finds **no** relation of any of $G$'s fold periods
$A_B,A_L,A_D,B_B,B_L,B_D,C_B,C_L,C_D$ with $\{m(2)=L'(E,0),\ \Omega_{12},\ \eta_{12},\ \Omega_{23},\ \eta_{23},\ \pi,\ 1\}$
at `maxcoeff` $10^6$, tolerance $10^{-70}$. So the Mahler measure of $G$'s curve
does **not** enter the fold periods linearly, despite the isogeny.

---

## 6. Hadamard finite part $\mathrm{FP}\int_0^1 H/(1-t)\,dt$  (Task 7)

The path is $[0,1)$ with $H$ continued from the upper half plane (§0), so the value is complex.

**For $E_m,F_m$** the puncture is a branch point: $H(t)=g(t)(1-t)^{-1/2}$ near $1$ with
$g(1)=-1/\sqrt{(a-1)(b-1)}=-1/\sqrt{\Delta}$, and
$$\mathrm{FP}:=\lim_{x\to1^-}\Bigl[\int_0^x\frac{H\,dt}{1-t}\;-\;2g(1)\,(1-x)^{-1/2}\Bigr].$$
**For $G$** the puncture is regular, $H(1)=-i/\sqrt{231}$, and the divergence is logarithmic:
$$\mathrm{FP}:=\lim_{x\to1^-}\Bigl[\int_0^x\frac{H\,dt}{1-t}\;+\;H(1)\log(1-x)\Bigr].$$
(For $G$ this limit was independently confirmed by direct quadrature at
$x=0.999,\,0.99999,\,0.9999999999$, converging to the stated value.)

| host | $g(1)$ or $H(1)$ | $\operatorname{Re}\mathrm{FP}$ | $\operatorname{Im}\mathrm{FP}$ |
|---|---|---|---|
| $E_1$ | $-0.20412414523193150818$ | −0.2730007970606116364225303671494942746988 | 0.7106476531948004379527235417391954753394 |
| $E_2$ | $-0.13363062095621219234$ | −0.2253096519629585879104655312560044678743 | 0.4470406597893019698378782311893503853748 |
| $E_3$ | $-0.10660035817780521715$ | −0.1985860641750073629548995073563698659261 | 0.3525032159685346419719436159465992345884 |
| $E_4$ | $-0.091287092917527685576$ | −0.1805110192914030513562807865814979489478 | 0.3002010941439971990698319635129378782194 |
| $E_5$ | $-0.081110710565381269082$ | −0.1671020422092576434395996418921646185715 | 0.2658736499273602635622045542444299124309 |
| $E_6$ | $-0.073720978077448566729$ | −0.1565823139280965490520887933321900783078 | 0.2411377050297526594010695344451286357159 |
| $E_{12}$ | $-0.051571062312939670362$ | −0.1205008176650495899575436492842283484563 | 0.1678109763415528376784260608001629847455 |
| $F_1$ | $-0.11785113019775792073$ | −0.1908404847834668001069148167422075415336 | 0.4030782993877955175468247848286412908321 |
| $F_2$ | $-0.077151674981045955131$ | −0.1595242417251179309723147341493185927980 | 0.2532796811970979500227194357376116495045 |
| $F_3$ | $-0.061545745489666366367$ | −0.1416996785139270462414363698303150371554 | 0.1996506236671019192794550948718790272813 |
| $G$ | $-i\cdot0.065795169495976898838$ | −0.3696835269269747340586558070193370117004 | 0.0516712187014800771526474971103900284646 |

**Every finite part is nonzero.**

**Requested PSLQ against $\{1,\Omega,\pi,\log2,\log3,\log m\}$:** **no relation** for any host
(`maxcoeff` $10^6$, tolerance $10^{-70}$). The only hits returned were degenerate relations
*inside the basis* (e.g. $\log2-\log m=0$ when $m=2$), with coefficient $0$ on the target.
**[numerical, 70 digits]**

**But the finite part is completely identified by the second-kind reduction (IV):**

> **(VI)** For $E_m,F_m$:
> $$\mathrm{FP}=f\Bigl[(\Omega_{01}-\eta_{01})-(\Omega_{23}-\eta_{23})\Bigr]-\frac{2}{\Delta}\;+\;i\,f\,(\Omega-\eta).$$

**[verified numerically, all ten hosts, residual $\lesssim10^{-90}$, `fpscan.out`]**.
In particular $\operatorname{Im}\mathrm{FP}=B_B=|\Delta_B|$ exactly — the finite part of the
integral to the puncture *sees the fold-difference period as its imaginary part*.
Note that $\Omega_{23}=\omega_1(E)$ is the real period, so (VI) expresses $\mathrm{FP}$
entirely in periods and quasi-periods of the curve plus the rational $-2/\Delta$.

For $G$ no analogous elementary identification was found.

---

## 7. What the CDT ($\tau=3/2$) theorem says in closed form (Task 6)

Take $E_m$ or $F_m$. The two folds give
$$c_{\rm near}[k]=A[k]\ (\text{real}),\qquad c_{\rm far}[k]=A[k]+i\,B[k].$$
A rational vector $(\alpha,\beta,\delta)$ with $\alpha+\beta c_i[B]+\delta c_i[D]=0$ for **both**
folds is, separating real and imaginary parts, exactly a solution of the **pair**
$$\beta\,B_B+\delta\,B_D=0,\qquad \alpha+\beta\,A_B+\delta\,A_D=0 .$$
Substituting (I) and (III), the first equation reads
$$\Bigl[\beta f+\delta\frac{a+b}{\Delta}\Bigr]\Omega\;-\;\delta\frac{2ab}{\Delta}\,\eta\;+\;\frac{\delta f}{2}\,(\Omega-\eta)\log c\;=\;0 .$$

* If $\delta=0$ then $\beta f(\Omega-\eta)=0$, impossible since $\Omega>\eta>0$; so $\beta=0$, then $\alpha=0$.
* If $\delta\ne0$, dividing by $\delta f/2$ (and using $f=ab/\Delta$, so that the $\eta$-coefficient
  becomes exactly $-4$) gives
  $$\boxed{\ \log c\;=\;\frac{p\,\Omega+4\,\eta}{\Omega-\eta}\ },\qquad
    p=-\Bigl(\frac{2\beta}{\delta}+\frac{2(a+b)}{ab}\Bigr)\in\mathbb{Q}\ \text{arbitrary}.$$

So **the CDT $\tau=3/2$ statement for these hosts is precisely**:

> the number $\log\dfrac{(a-1)(b-1)}{ab}$ is not equal to $\dfrac{p\,\Omega+4\,\eta}{\Omega-\eta}$
> for any rational $p$ (together with the accompanying rationality of $\alpha$),
> where $\Omega=\int_{r_1}^{r_2}\!\frac{dx}{y}$ and $\eta=\int_{r_1}^{r_2}\!\frac{x\,dx}{y}$
> are the period and quasi-period of $E$ over the fold cycle.
> The value of $p$ that would be required is $p_0=\bigl[(\Omega-\eta)\log c-4\eta\bigr]/\Omega$.
> **[numerical, 25 digits]** None is a small-height rational:
>
> | host | $p_0$ | host | $p_0$ |
> |---|---|---|---|
> | $E_1$ | −1.059777646043084742759739 | $E_6$ | −0.454903328211256307294003 |
> | $E_2$ | −0.693918832101824780716914 | $E_{12}$ | −0.395673087552090762372747 |
> | $E_3$ | −0.573975482477393022132573 | $F_1$ | −0.872757338558432040264136 |
> | $E_4$ | −0.514334807300505787272097 | $F_2$ | −0.491853224284370347783247 |
> | $E_5$ | −0.478651336473839841337950 | $F_3$ | −0.367101360313311144975849 |
>
> PSLQ of $p_0$ against $1$ with `maxcoeff` $=10^{12}$ and tolerance $10^{-70}$ returns **None**
> for all ten hosts: no rational of height $\le10^{12}$ works. **[numerical, 70 digits]**

Equivalently: **$\{\Omega,\ \eta,\ \Omega\log c,\ \eta\log c\}$ admits no $\mathbb{Q}$-linear
relation of the two-parameter shape above.** Concretely the ten logarithms involved are
$$\log\tfrac23,\ \log\tfrac79,\ \log\tfrac{22}{27},\ \log\tfrac56,\ \log\tfrac{38}{45},\ \log\tfrac{23}{27},\ \log\tfrac{47}{54}\ (E_m);\qquad
\log\tfrac{18}{25},\ \log\tfrac{21}{25},\ \log\tfrac{22}{25}\ (F_m).$$

For $G$ the same reduction does **not** close: there are three folds, the differences are
$c_2-c_1=i B$ and $c_3-c_2=-C$, and $B_B,C_B,B_D,C_D$ satisfy no $\mathbb{Q}$-linear relation with
$\Omega_{ij},\eta_{ij},\pi,L'(E,0)$; only the bilinear (V) is available. Hence for $G$ the
linear-independence statement remains genuinely two-dimensional and is **not** reducible to
a single logarithm.

### Summary table

| host | $\Delta_B$ identified? | $\Delta_D$ identified? | $\Delta_L$ identified? |
|---|---|---|---|
| $E_1,\dots,E_{12}$, $F_1,F_2,F_3$ | **yes** — $\pm i f(\Omega-\eta)$ | **yes** — $\pm i\bigl[\frac{(a+b)\Omega-2ab\eta}{\Delta}+\frac f2(\Omega-\eta)\log c\bigr]$ | **yes** — $\pm i\frac{\Omega}{2}\log c$ |
| $G$ | only via the bilinear (V) | no | no |
