<!-- LIAISON WORKFLOW (for any Claude session acting as certificate liaison):
River downloads files from their ChatGPT account (GPT Sol's documents, certificates, scripts) and tells you the path.
For each file: inspect; copy into certificates/<strand>/ (strands: book, eisenstein, domb, harmonic, catalan, bz-zeta5, padic, zeta5-zeta7, misc);
add a README line (what it is, which ledger ID/claim it supports); run anything runnable with `timeout 110` (python3/gp/lean);
log PASS/FAIL or output summary in certificates/INDEX.md; update the status column below. Never delete; never overwrite repo files;
do not touch consolidation/ or lattice/ (Fable's). Keep replies short: received / filed / ran / still wanted. -->

# Wanted list — artefacts to pull from ChatGPT ("GPT Sol")

Status legend: WANTED / RECEIVED / FILED / NOT NEEDED (see certificates/INDEX.md for filed items).

## P0 — highest priority

1. **The 15 canonical TeX sources** (ledger Source Index sheet) + any accompanying scripts/JSON certificates:
   - main(10).tex — Modular anchors for Apéry-like companions
   - main(20260813-131951).tex — Eisenstein sources of the sporadic companions
   - main(6).tex — A manufactured Apéry apparatus for a cusp-form L-value
   - main(8).tex — Harmonic Jets of Hypergeometric Sums
   - main(4).tex — Harmonic companion rows for binomial powers
   - main(7).tex — A two-level digit law for Apéry's pair
   - main(20260813-132000).tex — The fifteen sporadic Apéry-like pairs
   - main(3).tex — Why Apéry is alone
   - CATALAN_LITERAL_WINDOW_PAPER.tex — literal-window construction-quality bound (0.581983…)
   - CATALAN_LITERAL_WINDOW_PROOF.tex — long-form proof companion
   - frobenius.tex — Frobenius constants / middle-root obstruction
   - padic-apery-limits.tex — p-adic Apéry limits programme
   - main(9).tex — Foundations of the Λ-algebra
   - main(5).tex — The companion problem, narrative/provenance survey
   **Status: RECEIVED + FILED 2026-08-21 — all 14 listed items, all compile clean. See INDEX.md.**
   (Header says 15; only 14 are enumerated. Confirm whether one source is missing from the list.)

2. **Catalan "1-ε" paper/notes + ongoing Lean4 formalization** (mentioned yesterday, unsaved chats). Status: WANTED

3. **Bogner (n+1)|A_n proof notes** (yesterday). Status: WANTED

4. **Eisenstein Source Theorem**: Sturm-bound certificate scripts for the twelve sporadic source identities (ledger E01/E02), and the Lean project.
   Note: a Lean tarball already exists **mis-filed** under `packages/phase 0 math campaign/catalan/catalan two row 4e730650-0bf2-4257-baaf-25637b51474a-aristotle.tar(.gz)` — do not move it, just flag to Fable. Status: WANTED (scripts), NOTED (mis-filed tarball found in repo, not touched)

5. **Domb cuspidal apparatus certificates**: E07 (q^40 identity), E09 (denominators), E10 (L(f6,3)/2 PSLQ) — scripts. Status: WANTED

6. **β(4) level-24 dossier**: degree-4 elimination certificate, rank-16 cyclic-vector certificate mod 1000003 (cf. book/v8/08b). Status: WANTED

7. **L(3,χ5) dossier exact certificates**: order-6 Hadamard operator, branch matrix. Status: WANTED

8. **Level-100 Franke–η cusp-constant matrix computation**. Status: WANTED

## P1

- Shvets compression / L(4,χ-3) level-24 pair scripts
- "RM Catalan–Zeta Atlas" and "Sporadic Apéry Atlas" documents
- ζ(5) X0+(169) finite-field no-go certificates
- Γ1(7) Plücker computations
- Level-144 Catalan dossier scripts
- BZ κ-vector 434-digit PSLQ scripts
- s7 Lucas completion scripts
- Harmonic-jets Ore-algebra certificates (H05)
- ~~Modular parametrisation of the (12,4,16) third-order row (Almkvist–Zudilin "T")~~ — **NO LONGER NEEDED: derived independently 2026-08-21**, level 8, $t=(\eta_1\eta_8/\eta_2\eta_4)^8$, $F=-\tfrac16E_2(\tau)+\tfrac16E_2(2\tau)-\tfrac13E_2(4\tau)+\tfrac43E_2(8\tau)$ on $\Gamma_0(8)$. See `consolidation/ZETA3_TWO_LATTICE.md` §9. Still worth collecting GPT Sol's notes for cross-check.

All: Status WANTED (nothing received yet as of 2026-08-21).
