# nome_integrality — scripts for `consolidation/NOME_INTEGRALITY.md`

| script | checks |
|---|---|
| `sig.gp`      | Lemma A: `sqrt(1-34t+t^2)` integral; `v_2` profile |
| `basic.gp`    | Apéry recurrence vs binomial sum; `F*sigma`, `q`, `t(q)`, `F(t(q))` integral to `t^120` |
| `dwork2.gp`   | log-solution `G` from the ODE; identity `theta(G/F) = 1/(F sigma) - 1` to `t^400` |
| `dwork3.gp`   | (D1) `F(t)/F(t^p)` in `Z_p[[t]]`; (D2) `u(t)-u(t^p)/p` in `t Z_p[[t]]`; Dwork's `h` in `1+ptZ_p[[t]]` |
| `gauss.gp`    | Theorem N(2): Gauss congruences for `r_n` (`p<=60`, `n<=1200`) |
| `nome_big.gp` | same at `N=4000`, `p<=200`  (run on snake) |
| `rev.gp`      | `t(q)`, `F(t(q))`, `F*t` integral to `q^1200`  (run on snake) |
| `dwtrunc.gp`  | Theorem 8.3 truncated Dwork congruences, `p=2,3,5,7`, `n<=2600` |
| `kr.gp`       | Lemma 7.2: `B_n` (ODE) = `2*sum_j C(n,j)^2 C(n+j,j)^2 (H_{n+j}-H_{n-j})`, `n<=60`; Lemma 7.1 |
| `cn.gp`       | Theorem N(3): the exponents `c_n` of `q/t = prod (1-t^n)^{-c_n}` |
| `ct2.py`      | Theorem 8.1 `A_n = CT(Lambda^n)`, `n<=30`; Fact 8.2 reflexivity of the Newton polytope |
