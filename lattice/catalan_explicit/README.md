# lattice/catalan_explicit — explicit moment constructions

Write-up: `consolidation/CATALAN_EXPLICIT.md`.

`moments.gp` is the core: `nestgen(m,j)` returns `[A,B]` with
`int_0^1 int_0^1 K_Z^{(m)}(x,y) * (xy/(1-xy))^j dx dy = A*G - B`, exactly, in Q.
It generalises `nest()` from `lattice/catalan_positivity/rows_common.gp`
(the case `(m,j) = (3n,n)`); `j=0` is the Zudilin row at index m, `j=n` the Nesterenko (4,7) row.
Valid for `0 <= j <= m`; the integral diverges for `j > m`.

Run any script as
    cat ../catalan_positivity/rows_common.gp moments.gp X.gp > run.gp; gp -q run.gp
(`signs.gp` needs only `rows_common.gp`).

| script | what it measures |
|---|---|
| `momtab.gp`   | den(A),den(B), v2, sizes for all (n,j), n<=5 |
| `sos.gp`      | Hankel determinants, min over integer Q of Q^t H Q, lcm-denominator objective |
| `sossearch.gp`| box search over P = Q^2 and w*Q^2 with TRUE post-cancellation denominators |
| `signs.gp`    | sign pattern / 2-adic structure / ratio law of the cone-minimum coefficients, n<=44 |
| `control.gp`  | 2-row cone-min at small n; full (3n+1)-dim moment lattice at modulus 2^t |
| `pairs.gp`    | canonical (modulus-free) cone minima over ALL moment pairs |
| `bridge.gp`   | v2 of the mixed minor and the half-log-covolume predictor F over all pairs |
