import csv

from pristine.panel import diagnostic_sites, merge_panels, orient, panel_records, write_diagnostic_sites
from pristine.seqio import read_fasta, revcomp, write_fasta

REF = "ACGTTGCAGTCCATGGAAAAACTGACGTTAGCCGTACGATCGGATCCTAGCTAGGCTTACGG"


def _mut(seq, pos, base):
    return seq[:pos] + base + seq[pos + 1:]


def _site(sites, position, kind):
    return next(s for s in sites if s["position"] == position and s["type"] == kind)


def test_panel_one_record_per_distinct_sequence_and_species():
    other = _mut(REF, 10, "A")
    rows = [("Sp_a_1", "Sp_a", "target", REF), ("Sp_a_2", "Sp_a", "target", REF.lower()),
            ("bg_1", "bg", "background", other), ("Sp_b_1", "Sp_b", "non-target", other)]
    recs = panel_records("loc1", rows)
    assert [h for h, _ in recs] == ["loc1|species=Sp_a|genome=Sp_a_1|role=target",
                                    "loc1|species=Sp_b|genome=Sp_b_1|role=non-target",
                                    "loc1|species=bg|genome=bg_1|role=background"]
    assert all(s == s.upper() for _, s in recs)


def test_orient_returns_reference_orientation():
    assert orient(REF, revcomp(_mut(REF, 5, "A"))) == _mut(REF, 5, "A")


def test_substitution_and_indels_with_homopolymer_context():
    nt = _mut(REF, 10, "A")                       # substitution at position 11
    nt = nt[:17] + nt[18:]                        # 1-nt deletion in the AAAAA run (ref 16-20, 0-based 15-19)
    bg = REF[:40] + REF[44:]                      # 4-nt deletion "GGAT" at positions 41-44
    bg2 = REF[:30] + "TTT" + REF[30:]             # insertion after position 30
    rows = [("T_1", "T", "target", REF), ("N_1", "N", "non-target", nt),
            ("B_1", "B", "background", bg), ("C_1", "C", "background", bg2)]
    species, sites = diagnostic_sites(REF, rows)
    assert species == ["T", "N", "B", "C"]
    sub = _site(sites, 11, "substitution")
    assert (sub["ref"], sub["T"], sub["N"], sub["separates_target_from"]) == ("C", "C", "A", "N")
    hp = next(s for s in sites if s["type"] == "indel" and s["N"] == "-" and s["ref"] != "-")
    assert (hp["length"], hp["homopolymer_base"], hp["homopolymer_length"], hp["weak"]) == (1, "A", 5, "yes")
    dele = next(s for s in sites if s["type"] == "indel" and s["B"] == "-" and s["ref"] != "-")
    assert (dele["length"], dele["homopolymer_length"], dele["weak"]) == (4, 1, "no")
    ins = next(s for s in sites if s["type"] == "indel" and s["ref"] == "-")
    assert ins["C"] == "TTT" and ins["T"] == "-" and ins["separates_target_from"] == "C"
    assert len(sites) == 4


def test_polymorphic_species_still_separates_when_no_state_is_shared():
    rows = [("T_1", "T", "target", REF), ("N_1", "N", "non-target", _mut(REF, 3, "A")),
            ("N_2", "N", "non-target", _mut(REF, 3, "G"))]
    _, sites = diagnostic_sites(REF, rows)
    assert len(sites) == 1 and sites[0]["N"] == "A/G" and sites[0]["separates_target_from"] == "N"


def test_polymorphic_within_species_is_not_diagnostic_alone(tmp_path):
    rows = [("T_1", "T", "target", REF), ("T_2", "T", "target", _mut(REF, 3, "A")),
            ("N_1", "N", "non-target", REF)]
    out = tmp_path / "s.tsv"
    assert write_diagnostic_sites(out, REF, rows) == 0
    assert next(csv.reader(open(out), delimiter="\t"))[:3] == ["position", "type", "length"]


def test_merge_prefixes_colliding_loci(tmp_path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    write_fasta(a, [("loc1|species=X|role=target", "ACGT"), ("loc2|species=X|role=target", "AC")])
    write_fasta(b, [("loc1|species=Y|role=target", "GGGG")])
    merge_panels({"ga": a, "gb": b, "gc": tmp_path / "missing.fa"}, tmp_path / "m.fa")
    heads = [h for h, _ in read_fasta(tmp_path / "m.fa")]
    assert heads == ["ga_loc1|species=X|role=target", "loc2|species=X|role=target", "gb_loc1|species=Y|role=target"]


def test_weak_homopolymer_indels_are_not_counted():
    from pristine.pipeline import _events
    from pristine.scoring import sites_between, _state_matrix, _or_rows

    hp = REF[:17] + REF[18:]                      # 1-nt deletion in AAAAA
    out = REF[:3] + REF[4:]                       # 1-nt deletion outside homopolymers (T in GTT? run 2)
    assert _events(REF, hp) == 0
    assert _events(REF, REF[:17] + "A" + REF[17:]) == 0   # 1-nt insertion into AAAAA
    assert _events(REF, out) == 1
    assert _events(REF, REF[:40] + REF[44:]) == 1          # multi-nt deletion counts

    # alignment scoring: column 17 of the AAAAA run gapped in the other group
    aln = [REF, REF[:17] + "-" + REF[18:], REF[:5] + "-" + REF[6:]]
    mat = _state_matrix(aln)
    s = sites_between(_or_rows(mat, [0]), _or_rows(mat, [1]), aln)
    assert s.count == 0 and s.weak_indels == [17]
    s = sites_between(_or_rows(mat, [0]), _or_rows(mat, [2]), aln)
    assert s.count == 1 and s.weak_indels == []
