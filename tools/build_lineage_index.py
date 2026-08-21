#!/usr/bin/env python3
"""Reconstruct branch lineage across the ChatGPT research archive.

Each conversation export's title looks like:
  "<project> - Branch · Branch · <base title>"
Depth = number of "Branch · " occurrences. Each depth level is a genuinely
separate ChatGPT conversation (distinct URL) with no shared message content
with its parent, so lineage must be inferred rather than read off the file.

Heuristic: group files by (project, base_title). Within a group, sort by
depth. For each node at depth d, its parent is the node at depth d-1 in the
same group with the largest capturedAt that is still <= this node's
capturedAt (i.e. the most recently-branched-from ancestor available at the
time this branch was made). If no such node exists, fall back to the
earliest depth-(d-1) node, and if still none, it's a root.

Output: a JSON index (family -> ordered chain per leaf) and a Markdown
summary, written to chatgpt-research-archive/lineage_index.{json,md}.
"""
import json
import re
import glob
import os
from collections import defaultdict

ARCHIVE_ROOT = "chatgpt-research-archive/conversations"
OUT_JSON = "chatgpt-research-archive/lineage_index.json"
OUT_MD = "chatgpt-research-archive/lineage_index.md"

BRANCH_RE = re.compile(r"^(.*?) - ((?:Branch · )*)(.*)$")


def parse_title(title):
    m = BRANCH_RE.match(title)
    if not m:
        return None, 0, title
    project, branch_part, base = m.groups()
    depth = branch_part.count("Branch · ")
    return project, depth, base


def load_records():
    records = []
    for path in glob.glob(f"{ARCHIVE_ROOT}/*/*.json"):
        with open(path) as f:
            d = json.load(f)
        title = d.get("title", "")
        project, depth, base = parse_title(title)
        records.append({
            "path": path,
            "project_dir": os.path.basename(os.path.dirname(path)),
            "project_title": project,
            "depth": depth,
            "base_title": base,
            "capturedAt": d.get("capturedAt"),
            "url": d.get("url"),
            "n_messages": len(d.get("messages", [])),
        })
    return records


def build_groups(records):
    groups = defaultdict(list)
    for r in records:
        key = (r["project_dir"], r["base_title"])
        groups[key].append(r)
    for key in groups:
        groups[key].sort(key=lambda r: (r["depth"], r["capturedAt"] or ""))
    return groups


def assign_parents(group):
    by_depth = defaultdict(list)
    for r in group:
        by_depth[r["depth"]].append(r)
    for r in group:
        r["parent_path"] = None
        if r["depth"] == 0:
            continue
        candidates = by_depth[r["depth"] - 1]
        if not candidates:
            continue
        before = [c for c in candidates if (c["capturedAt"] or "") <= (r["capturedAt"] or "")]
        parent = max(before, key=lambda c: c["capturedAt"] or "") if before else min(
            candidates, key=lambda c: c["capturedAt"] or ""
        )
        r["parent_path"] = parent["path"]
    return group


def build_chains(group):
    """For each leaf (node with no children), return its root->leaf chain of paths."""
    by_path = {r["path"]: r for r in group}
    children = defaultdict(list)
    for r in group:
        if r["parent_path"]:
            children[r["parent_path"]].append(r["path"])
    leaves = [r["path"] for r in group if r["path"] not in children]
    chains = []
    for leaf in leaves:
        chain = []
        cur = leaf
        seen = set()
        while cur and cur not in seen:
            seen.add(cur)
            chain.append(cur)
            cur = by_path[cur]["parent_path"]
        chain.reverse()
        chains.append(chain)
    return chains


def main():
    records = load_records()
    groups = build_groups(records)

    index = {}
    md_lines = ["# Lineage index\n", "Reconstructed branch families across the archive.\n"]

    for (project_dir, base_title), group in sorted(groups.items()):
        group = assign_parents(group)
        chains = build_chains(group)
        family_key = f"{project_dir}::{base_title}"
        index[family_key] = {
            "project_dir": project_dir,
            "base_title": base_title,
            "nodes": [
                {
                    "path": r["path"],
                    "depth": r["depth"],
                    "capturedAt": r["capturedAt"],
                    "parent_path": r["parent_path"],
                    "n_messages": r["n_messages"],
                }
                for r in group
            ],
            "chains": chains,
        }

        md_lines.append(f"\n## {project_dir} / {base_title}\n")
        md_lines.append(f"- {len(group)} node(s), {len(chains)} leaf chain(s)\n")
        for chain in chains:
            names = [os.path.basename(p) for p in chain]
            md_lines.append(f"  - " + " -> ".join(names))

    with open(OUT_JSON, "w") as f:
        json.dump(index, f, indent=2)
    with open(OUT_MD, "w") as f:
        f.write("\n".join(md_lines) + "\n")

    total_families = len(index)
    total_files = sum(len(v["nodes"]) for v in index.values())
    total_chains = sum(len(v["chains"]) for v in index.values())
    print(f"Families: {total_families}, files: {total_files}, leaf chains: {total_chains}")
    print(f"Wrote {OUT_JSON} and {OUT_MD}")


if __name__ == "__main__":
    main()
