import shutil
from io import StringIO
from pathlib import Path

import pytest

from pristine import blast

XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastn</BlastOutput_program>
  <BlastOutput_version>BLASTN 2.16.0+</BlastOutput_version>
  <BlastOutput_reference>x</BlastOutput_reference>
  <BlastOutput_db>nt</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>locusA</BlastOutput_query-def>
  <BlastOutput_query-len>1000</BlastOutput_query-len>
  <BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect></Parameters></BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-def>locusA</Iteration_query-def>
      <Iteration_query-len>1000</Iteration_query-len>
      <Iteration_hits>
        <Hit><Hit_num>1</Hit_num><Hit_id>gi|1|gb|CP000001.1|</Hit_id>
          <Hit_def>Yersinia pseudotuberculosis strain Z chromosome, complete genome</Hit_def>
          <Hit_accession>CP000001</Hit_accession><Hit_len>4700000</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>1800</Hsp_bit-score><Hsp_score>990</Hsp_score>
            <Hsp_evalue>0</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>1000</Hsp_query-to>
            <Hsp_hit-from>10</Hsp_hit-from><Hsp_hit-to>1009</Hsp_hit-to><Hsp_query-frame>1</Hsp_query-frame>
            <Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>1000</Hsp_identity><Hsp_positive>1000</Hsp_positive>
            <Hsp_gaps>0</Hsp_gaps><Hsp_align-len>1000</Hsp_align-len><Hsp_qseq>A</Hsp_qseq><Hsp_hseq>A</Hsp_hseq><Hsp_midline>|</Hsp_midline></Hsp></Hit_hsps></Hit>
        <Hit><Hit_num>2</Hit_num><Hit_id>gi|2|gb|CP000002.1|</Hit_id>
          <Hit_def>Yersinia similis strain Q chromosome</Hit_def>
          <Hit_accession>CP000002</Hit_accession><Hit_len>4700000</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>1700</Hsp_bit-score><Hsp_score>900</Hsp_score>
            <Hsp_evalue>0</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>1000</Hsp_query-to>
            <Hsp_hit-from>10</Hsp_hit-from><Hsp_hit-to>1009</Hsp_hit-to><Hsp_query-frame>1</Hsp_query-frame>
            <Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>985</Hsp_identity><Hsp_positive>985</Hsp_positive>
            <Hsp_gaps>0</Hsp_gaps><Hsp_align-len>1000</Hsp_align-len><Hsp_qseq>A</Hsp_qseq><Hsp_hseq>A</Hsp_hseq><Hsp_midline>|</Hsp_midline></Hsp></Hit_hsps></Hit>
        <Hit><Hit_num>3</Hit_num><Hit_id>gi|3|gb|CP000003.1|</Hit_id>
          <Hit_def>Yersinia pestis CO92 chromosome</Hit_def>
          <Hit_accession>CP000003</Hit_accession><Hit_len>4600000</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>1600</Hsp_bit-score><Hsp_score>850</Hsp_score>
            <Hsp_evalue>0</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>1000</Hsp_query-to>
            <Hsp_hit-from>10</Hsp_hit-from><Hsp_hit-to>1009</Hsp_hit-to><Hsp_query-frame>1</Hsp_query-frame>
            <Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>960</Hsp_identity><Hsp_positive>960</Hsp_positive>
            <Hsp_gaps>0</Hsp_gaps><Hsp_align-len>1000</Hsp_align-len><Hsp_qseq>A</Hsp_qseq><Hsp_hseq>A</Hsp_hseq><Hsp_midline>|</Hsp_midline></Hsp></Hit_hsps></Hit>
        <Hit><Hit_num>4</Hit_num><Hit_id>gi|4|gb|CP000004.1|</Hit_id>
          <Hit_def>Escherichia coli K-12 partial sequence</Hit_def>
          <Hit_accession>CP000004</Hit_accession><Hit_len>4600000</Hit_len>
          <Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>300</Hsp_bit-score><Hsp_score>150</Hsp_score>
            <Hsp_evalue>1e-50</Hsp_evalue><Hsp_query-from>1</Hsp_query-from><Hsp_query-to>300</Hsp_query-to>
            <Hsp_hit-from>10</Hsp_hit-from><Hsp_hit-to>309</Hsp_hit-to><Hsp_query-frame>1</Hsp_query-frame>
            <Hsp_hit-frame>1</Hsp_hit-frame><Hsp_identity>297</Hsp_identity><Hsp_positive>297</Hsp_positive>
            <Hsp_gaps>0</Hsp_gaps><Hsp_align-len>300</Hsp_align-len><Hsp_qseq>A</Hsp_qseq><Hsp_hseq>A</Hsp_hseq><Hsp_midline>|</Hsp_midline></Hsp></Hit_hsps></Hit>
      </Iteration_hits>
      <Iteration_stat><Statistics><Statistics_db-num>1</Statistics_db-num><Statistics_db-len>1</Statistics_db-len>
        <Statistics_hsp-len>0</Statistics_hsp-len><Statistics_eff-space>1</Statistics_eff-space>
        <Statistics_kappa>0.4</Statistics_kappa><Statistics_lambda>0.6</Statistics_lambda><Statistics_entropy>0.8</Statistics_entropy></Statistics></Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def test_guess_organism():
    assert blast.guess_organism("Yersinia pestis strain CO92 chromosome") == "Yersinia pestis"
    assert blast.guess_organism("PREDICTED: Bos taurus mRNA") == "Bos taurus"
    assert blast.guess_organism("Yersinia sp. FDAARGOS_228 chromosome") == "Yersinia sp."
    assert blast.guess_organism("MAG: Candidatus Methanoplasma termitum genome") == "Methanoplasma termitum"
    assert blast.guess_organism("[Clostridium] innocuum strain X") == "Clostridium innocuum"
    assert blast.guess_organism("uncultured bacterium clone 7") == "uncultured"


def test_parse_and_summarize_blast_xml():
    hits = blast.parse_blast_xml(StringIO(XML), min_identity=90.0)
    assert list(hits) == ["locusA"]
    assert [h.organism for h in hits["locusA"]] == ["Yersinia pseudotuberculosis", "Yersinia similis", "Yersinia pestis", "Escherichia coli"]
    assert hits["locusA"][1].pct_identity == 98.5 and hits["locusA"][3].query_coverage == 30.0
    # 30 informative sites in a 1000 bp amplicon -> tolerated identity 97 %
    tol = {"locusA": blast.tolerated_identity(1000, 30, fallback=90.0)}
    assert tol["locusA"] == 97.0
    table = blast.summarize_blast(hits, ["Yersinia_pseudotuberculosis"], ["Yersinia_pestis"], tol)
    status = dict(zip(table.Organism, table.Status))
    assert status["Yersinia pseudotuberculosis"] == blast.STATUS_TARGET
    assert status["Yersinia pestis"] == blast.STATUS_NONTARGET
    assert status["Yersinia similis"] == blast.STATUS_UNKNOWN
    risk = dict(zip(table.Organism, table.Risk))
    assert risk == {"Yersinia pseudotuberculosis": False, "Yersinia similis": True, "Yersinia pestis": False, "Escherichia coli": False}


def test_tolerated_identity_fallback():
    assert blast.tolerated_identity(1000, None, 98.0) == 98.0
    assert blast.tolerated_identity(1000, 5, 98.0) == 99.5
    assert blast.tolerated_identity(1000, 100, 98.0) == 98.0  # never below the fallback


@pytest.mark.skipif(shutil.which("minimap2") is None, reason="minimap2 required")
def test_background_identity(dataset, tmp_path):
    from pristine.seqio import read_fasta_dict, write_fasta

    root, files = dataset
    target = read_fasta_dict(files["targets"][0])["chr"]
    amp = tmp_path / "amp.fasta"
    write_fasta(amp, [("locusA", target[20000:21200])])
    df = blast.background_identity(amp, [files["global_nt"][0], files["targets"][1]], tmp_path / "work", 1,
                                   {"locusA": 99.0})
    by = df.set_index("Background_Genome")
    assert 95 < by.loc["Alienus_tertius_1", "Pct_Identity"] < 99.0 and not by.loc["Alienus_tertius_1", "Risk"]
    assert by.loc["Targetus_primus_2", "Pct_Identity"] > 99.0 and by.loc["Targetus_primus_2", "Risk"]
