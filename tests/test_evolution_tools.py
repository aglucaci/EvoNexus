"""External-tool resolution for the evolutionary analyses.

The pipeline's honesty rule cuts both ways: a missing tool must be reported
as a documented exclusion, and a tool that is present must not be reported as
missing. Bioconda ships IQ-TREE 2 as `iqtree` in some builds and `iqtree2` in
others, which broke the second half of that rule.
"""

from __future__ import annotations

import pytest

from evonexus.evolution import pipeline, tools


@pytest.fixture
def no_tools(monkeypatch):
    """Start from an environment where nothing is installed."""
    installed: set[str] = set()
    monkeypatch.setattr(tools, "tool_available", lambda name: name in installed)
    return installed


def test_resolve_iqtree_prefers_version_2_binary(no_tools):
    no_tools.update({"iqtree", "iqtree2"})
    assert tools.resolve_iqtree() == "iqtree2"


def test_resolve_iqtree_accepts_the_unsuffixed_binary(no_tools):
    no_tools.add("iqtree")
    assert tools.resolve_iqtree() == "iqtree"


def test_resolve_iqtree_reports_nothing_when_absent(no_tools):
    assert tools.resolve_iqtree() is None


def test_run_iqtree_reports_a_missing_tool_rather_than_approximating(no_tools):
    result = tools.run_iqtree("aln.fa", "prefix")
    assert result["status"] == "tool_missing"
    assert result["treefile"] is None if "treefile" in result else True


@pytest.mark.parametrize("installed", [{"iqtree"}, {"iqtree2"}, {"iqtree", "iqtree2"}])
def test_installed_iqtree_is_not_reported_missing(monkeypatch, installed):
    """The regression: `iqtree` alone used to be reported as `iqtree2` missing."""
    monkeypatch.setattr(tools, "tool_available", lambda name: name in installed)

    missing = pipeline.missing_external_tools()

    assert "iqtree2" not in missing
    # MAFFT and HyPhy are genuinely absent here and must still be reported.
    assert missing == ["hyphy", "mafft"]


def test_absent_iqtree_is_still_reported(no_tools):
    no_tools.update({"mafft", "hyphy"})
    assert pipeline.missing_external_tools() == ["iqtree2"]


def test_a_complete_toolchain_reports_nothing_missing(no_tools):
    no_tools.update({"mafft", "hyphy", "iqtree"})
    assert pipeline.missing_external_tools() == []
