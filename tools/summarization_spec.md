# Summarization pass — spec for cheap-model run over the archive

Input: `chatgpt-research-archive/lineage_index.json` (built by `tools/build_lineage_index.py`).
Goal: one structured summary per **leaf chain** (not per file), tagging results and
triaging them against `packages/Odd_Zeta_Research_Master_Ledger.xlsx`.

## Why per-chain, not per-file

Branch files share **zero text** with their parent (confirmed by earlier survey) —
each file only contains messages from its own divergence point onward. Summarizing
a deep branch file alone is context-blind. `lineage_index.json`'s `chains` field
gives, per leaf, the ordered list of file paths from root to leaf — read and
summarize that whole chain as one unit.

Of the 57 families: 51 are single-node (1 file = 1 chain, trivial case). 6 are real
fanouts (see below) where multiple chains share a trunk — summarize the shared
trunk once per family, then only the diverging tail per chain, to avoid redundant
reading:

- `apery-systems::TODO: Ramanujan machine` (9 nodes, 7 chains — biggest)
- `apery-systems::lit review research fanout`
- `apery-systems-2::MASTER 11`
- `apery-systems-4::Zaggy Sporadic Sequences Analysis`
- `apery-systems-p-3::SMITH FROBENIUS CELIING!`
- `apery-systems-p-3::YOU CAN HELP!`

## Output schema (one JSON object per leaf chain)

```json
{
  "family": "apery-systems-p-3::SMITH FROBENIUS CELIING!",
  "chain": ["path/to/root.json", "path/to/branch.json", "..."],
  "leaf_path": "path/to/leaf.json",
  "title_human": "short human-readable label for this thread",
  "claims": [
    {
      "statement": "exact or near-exact quote of the mathematical claim — do not paraphrase formulas",
      "object": "which constant/sequence/construction this is about (e.g. 'Catalan level-144', 'zeta(5) level-12 CM-isogeny')",
      "status": "proved | conjectural | experimental/numerical | abandoned | superseded-by-later-message-in-thread",
      "source_file": "which file in the chain this came from",
      "notes": "anything about how confident the thread itself seemed, dead ends, handoffs requested"
    }
  ],
  "connections_to_literature": ["e.g. Zudilin Catalan permutation-group construction", "Brown-Zudilin zeta(5)"],
  "open_threads_or_handoffs": ["explicit TODOs, 'needs handoff from X', unresolved questions left at the end of the chain"],
  "in_ledger": "yes | no | partially — best guess after checking Odd_Zeta_Research_Master_Ledger.xlsx Ledger + Do Not Overstate sheets by keyword",
  "ledger_notes": "if in_ledger is yes/partially, which ledger ID(s) it corresponds to; if the ledger's Do Not Overstate sheet says this claim is superseded, flag that explicitly here",
  "flags": ["e.g. 'DEAD BRANCH per filename', 'title says FAILED/deadend', 'looks numerically wrong, double check'"]
}
```

Notes on fields:
- `status` should reflect what the **thread itself** claims, not the summarizer's
  own judgment of correctness — that's what `in_ledger`/`ledger_notes` is for.
- Respect filename signals already present in the archive: files named with
  `dead-branch-`, `failed-`, `deadend` should get `status: "abandoned"` by default
  unless the content clearly overturns that.
- Multiple `claims` per chain are expected and fine — don't force one claim per chain.

## Prompt template — single-node chain (51 of 57 families)

```
You are tagging one research conversation from a math archive (geometric "modular
sources" theory classifying Apéry-like sequences; a "multiple-lattice method" for
irrationality/worthiness exponents; connections to Zudilin's work on Catalan's
constant). Read the attached file (a full conversation transcript, .md preferred
over .json for readability). Extract every distinct mathematical claim made,
whether proved, conjectural, numerical/experimental, or abandoned. Quote formulas
exactly — do not paraphrase or simplify them. Note explicit handoffs/TODOs left at
the end. Then check packages/Odd_Zeta_Research_Master_Ledger.xlsx (Ledger sheet and
Do Not Overstate sheet — read via openpyxl) for whether each claim is already
tracked, and if the Do Not Overstate sheet flags it as superseded, say so explicitly
with the corrected value. Output strictly the JSON schema in
tools/summarization_spec.md (one object, single-file chain).
```

## Prompt template — fanout family (6 families: trunk + tail)

```
This family branches: multiple ChatGPT conversations diverge from shared trunk
messages, but branch files contain NO overlapping text with their parent (they only
capture messages from the divergence point onward) — so read the trunk file(s) and
each leaf's tail file(s) together, in the order given by the chain, to reconstruct
context. First, summarize the shared trunk once (claims established before the
first divergence). Then for each leaf chain, summarize only what's new/different in
its diverging tail, referencing the trunk summary rather than repeating it. Same
extraction rules and output schema as the single-node prompt (one JSON object per
leaf chain, but tail-only claims may reference "see trunk" for repeated context).
Cross-check against the ledger the same way.
```

## Execution notes

- Prefer the `.md` export of each file over `.json` for the actual read (same
  content, much more token-efficient — no HTML/metadata noise); `.json` is only
  needed for `capturedAt`/`url` bookkeeping, already extracted into
  `lineage_index.json`.
- This is a good fit for a cheap/fast model — extraction + a fixed lookup against
  one xlsx file, high volume, low reasoning depth per item.
- Suggested output location: one file per chain under
  `chatgpt-research-archive/summaries/<family-slug>__<leaf-basename>.json`, so
  results are diffable and don't require re-running everything if one summary needs
  a redo.
- After the pass: a cheap follow-up script can aggregate `in_ledger: "no"` claims
  across all summaries into a single triage queue — that's the actual worklist for
  "what's real research sitting outside the ledger."
