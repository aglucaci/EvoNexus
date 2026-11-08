"""Host phylogeny: Newick parsing and the covariance matrix used by model M5."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, List

import numpy as np


class Node:
    __slots__ = ("name", "length", "children", "parent")

    def __init__(self, name: str = "", length: float = 0.0):
        self.name = name
        self.length = length
        self.children: List["Node"] = []
        self.parent: "Node" | None = None

    def is_leaf(self) -> bool:
        return not self.children


def parse_newick(text: str) -> Node:
    """Parse a Newick string into a tree of :class:`Node`.

    Supports branch lengths, quoted labels and internal node labels. Raises
    ``ValueError`` on malformed input rather than returning a partial tree.
    """
    s = text.strip()
    if not s:
        raise ValueError("empty Newick string")
    if s.endswith(";"):
        s = s[:-1]

    pos = 0

    def parse_node() -> Node:
        nonlocal pos
        node = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                child = parse_node()
                child.parent = node
                node.children.append(child)
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"malformed Newick near position {pos}")
        # label
        start = pos
        if pos < len(s) and s[pos] in "'\"":
            quote = s[pos]
            pos += 1
            while pos < len(s) and s[pos] != quote:
                pos += 1
            node.name = s[start + 1: pos]
            pos += 1
        else:
            while pos < len(s) and s[pos] not in ",():;":
                pos += 1
            node.name = s[start:pos].strip()
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and (s[pos].isdigit() or s[pos] in ".eE+-"):
                pos += 1
            try:
                node.length = float(s[start:pos])
            except ValueError:
                node.length = 0.0
        return node

    root = parse_node()
    if pos != len(s):
        raise ValueError(f"unparsed Newick remainder at position {pos}: {s[pos:pos+20]!r}")
    return root


def _leaves(node: Node) -> List[Node]:
    if node.is_leaf():
        return [node]
    out: List[Node] = []
    for child in node.children:
        out.extend(_leaves(child))
    return out


def _root_distance(node: Node) -> float:
    d = 0.0
    cursor = node
    while cursor.parent is not None:
        d += cursor.length
        cursor = cursor.parent
    return d


def _path_to_root(node: Node) -> List[Node]:
    path = [node]
    while path[-1].parent is not None:
        path.append(path[-1].parent)
    return path


def newick_to_covariance(
    newick: str,
    species: List[str] | None = None,
    scale_to_correlation: bool = True,
) -> Dict[str, Any]:
    """Brownian-motion covariance: shared root-to-MRCA path length.

    Species names are normalized by replacing underscores with spaces, which is
    how bat species usually differ between a tree file and a metadata table.
    """
    root = parse_newick(newick)
    leaves = _leaves(root)
    names = [leaf.name.replace("_", " ").strip() for leaf in leaves]
    lookup = {name: leaf for name, leaf in zip(names, leaves)}

    if species:
        wanted = [s for s in species if s in lookup]
        missing = [s for s in species if s not in lookup]
    else:
        wanted, missing = names, []

    n = len(wanted)
    C = np.zeros((n, n))
    D = np.zeros((n, n))
    paths = {s: {id(node): _root_distance(node) for node in _path_to_root(lookup[s])} for s in wanted}

    for i, si in enumerate(wanted):
        for j, sj in enumerate(wanted):
            shared = set(paths[si]) & set(paths[sj])
            mrca_depth = max((paths[si][k] for k in shared), default=0.0)
            C[i, j] = mrca_depth
            D[i, j] = (_root_distance(lookup[si]) - mrca_depth) + (
                _root_distance(lookup[sj]) - mrca_depth
            )

    if scale_to_correlation and n:
        d = np.sqrt(np.clip(np.diag(C), 1e-12, None))
        C = C / np.outer(d, d)
        np.fill_diagonal(C, 1.0)

    return {
        "species": wanted,
        "matrix": C.tolist(),
        "distance_matrix": D.tolist(),
        "missing_species": missing,
        "scaled_to_correlation": scale_to_correlation,
        "n_species": n,
        "note": (
            "Brownian-motion covariance from shared root-to-MRCA path length. "
            "Species absent from the tree are listed in `missing_species` and "
            "are modelled with exchangeable host effects instead."
        ),
    }


def covariance_from_file(path: str | Path, species: List[str] | None = None,
                         scale_to_correlation: bool = True) -> Dict[str, Any]:
    text = Path(path).read_text(encoding="utf-8")
    return newick_to_covariance(text, species, scale_to_correlation)
