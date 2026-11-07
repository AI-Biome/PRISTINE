import pytest

from pristine.consensus import iupac_consensus, masked_template
from pristine.manifest import Genome
from pristine.scoring import is_informative, score_alignment


def _g(gid, species, target):
    return Genome(gid, None, species, "target" if target else "group_nontarget", "G")


def test_is_informative():
    assert is_informative({"A"}, {"G"})
    assert is_informative({"A", "G"}, {"C", "T"})
    assert not is_informative({"A", "G"}, {"G"})
    assert not is_informative(set(), {"G"})
    assert is_informative({"A"}, {"-"})  # indel counts as a state


def test_score_alignment_external_and_internal():
    genomes = {
        "t1": _g("t1", "Sp_a", True),
        "t2": _g("t2", "Sp_a", True),
        "t3": _g("t3", "Sp_b", True),
        "n1": _g("n1", "Nt_x", False),
        "n2": _g("n2", "Nt_x", False),
        "m1": _g("m1", "Nt_y", False),
    }
    #            0123456789
    aln = {
        "t1": "ACGTACGTAC",
        "t2": "ACGTACGTAC",
        "t3": "ACGTACGTTC",  # pos 8: Sp_b differs from Sp_a -> internal site
        "n1": "ACGAACGTAC",  # pos 3: Nt_x differs from group
        "n2": "ACGAACGTAC",
        "m1": "ACGTACG-AC",  # pos 7: gap in Nt_y -> informative
    }
    sc = score_alignment("L", aln, genomes)
    assert sc.external_positions["Nt_x"] == [3]
    assert sc.external["Nt_x"].substitutions == [3] and sc.external["Nt_x"].indel_events == []
    assert sc.external_positions["Nt_y"] == [7]
    assert sc.external["Nt_y"].indel_events == [7] and sc.external["Nt_y"].substitutions == []
    # pos 8: group states {A, T} vs Nt_x {A} -> not disjoint -> not informative
    assert sc.internal_positions[("Sp_a", "Sp_b")] == [8]
    assert sc.internal_positions[("Sp_b", "Sp_a")] == [8]
    assert sc.aln_len == 10
    assert sc.min_external_count == 1
    rec = sc.to_record()
    assert rec["SNP_Pos_External_Nt_x"] == "3"
    assert rec["Subst_External_Nt_x"] == 1 and rec["Indel_External_Nt_y"] == 1
    assert rec["Abs_Internal_Sp_a_vs_Sp_b"] == 1
    assert sc.count_in_region(0, 5) == {"External_Nt_x": 1, "External_Nt_y": 0,
                                        "Internal_Sp_a_vs_Sp_b": 0, "Internal_Sp_b_vs_Sp_a": 0}
    assert 0 < sc.combined_score < 1


def test_score_alignment_rejects_ragged():
    genomes = {"t1": _g("t1", "A", True), "n1": _g("n1", "B", False)}
    with pytest.raises(ValueError):
        score_alignment("L", {"t1": "ACGT", "n1": "ACG"}, genomes)


def test_masked_template_and_consensus():
    seqs = ["ACGTAC-TA", "ACGTAC-TA", "ACGAACGTA"]
    assert masked_template(seqs) == "ACGNACNTA"
    assert iupac_consensus(seqs) == "ACGWACGTA"
    assert masked_template(["ACNT", "ACNT"]) == "ACNT"


def test_indel_run_counts_as_one_event():
    genomes = {"t1": _g("t1", "A", True), "n1": _g("n1", "B", False)}
    #              0123456789012
    aln = {"t1": "ACGT------ACG",
           "n1": "ACGTTTTTTTACC"}  # 6-column insertion in n1 + substitution at 12
    sc = score_alignment("L", aln, genomes)
    assert sc.external["B"].indel_events == [4]
    assert sc.external["B"].substitutions == [12]
    assert sc.external_positions["B"] == [4, 12]
    assert sc.external["B"].count == 2
    assert sc.to_record()["Abs_External_B"] == 2
    # two separate gap runs are two events
    aln = {"t1": "AC--GT--AC", "n1": "ACGGGTTTAC"}
    sc = score_alignment("L", aln, genomes)
    assert sc.external["B"].indel_events == [2, 6]
