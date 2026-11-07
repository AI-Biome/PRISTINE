import pytest

from pristine.pseudoalign import build_alignment, flip_cs, hit_states, parse_cs
from pristine.seqio import revcomp

WINDOW = "ACGTACGTAC"  # 10 bp


def test_parse_cs():
    assert parse_cs(":4*ag+cc-t:3") == [(":", "4"), ("*", "ag"), ("+", "cc"), ("-", "t"), (":", "3")]
    with pytest.raises(ValueError):
        parse_cs(":4xx")


def test_hit_states_plus_strand():
    # genome: positions 0-3 match, pos 4 window A -> genome G, window 5-6 (CG) missing in genome,
    # genome has extra TT before window pos 7, then 7-9 match
    states, ins = hit_states(WINDOW, 0, 10, "+", ":4*ga+cg-tt:3")
    assert "".join(states) == "ACGTG--TAC"
    assert ins == {7: "TT"}


def test_hit_states_partial_alignment_is_N():
    states, ins = hit_states(WINDOW, 2, 8, "+", ":6")
    assert "".join(states) == "NNGTACGTNN"


def test_hit_states_minus_strand_equals_plus_strand_result():
    # Build the same alignment seen from the reverse strand and check we recover it.
    plus_cs = ":4*ga+cg-tt:3"
    exp_states, exp_ins = hit_states(WINDOW, 0, 10, "+", plus_cs)
    # minus-strand cs: reversed ops with complemented bases
    ops = parse_cs(plus_cs)
    minus_cs = "".join(op + payload for op, payload in flip_cs(ops))
    assert minus_cs == ":3-aa+cg*ct:4"
    states, ins = hit_states(WINDOW, 0, 10, "-", minus_cs)
    assert states == exp_states and ins == exp_ins


def test_hit_states_length_mismatch_raises():
    with pytest.raises(ValueError):
        hit_states(WINDOW, 0, 10, "+", ":9")


def test_build_alignment_inserts_columns():
    aln = build_alignment(
        WINDOW,
        {
            "ref": (0, 10, "+", ":10"),
            "g1": (0, 10, "+", ":4*ga+cg-tt:3"),
            "g2": (0, 10, "+", ":7-t:3"),
        },
    )
    assert aln["ref"] == "ACGTACG--TAC"
    assert aln["g1"] == "ACGTG--TTTAC"
    assert aln["g2"] == "ACGTACGT-TAC"
    assert len({len(s) for s in aln.values()}) == 1
