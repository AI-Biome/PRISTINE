# PRISTINE: PRimer-based Identification Suite Targeting Informative Nucleotide Elements

PRISTINE designs PCR primer pairs for amplicon-sequencing based identification of a
target group of bacteria (one or several species) against closely related non-target
species. Each primer pair

- binds sites that are **identical in every input genome** (targets and non-targets),
  so the amplicon is produced from all of them;
- encloses a region that is **rich in informative sites** (positions where the target
  group differs from each non-target species), so reads of the amplicon tell targets
  and non-targets apart;
- yields **exactly one product per genome in silico**, and none in optional
  background genomes.

Amplicon sizes of 500-1800 bp are intended for Oxford Nanopore sequencing and for
downstream analysis with [PRISCREEN](https://github.com/AI-Biome/PRISCREEN).

## Installation

Dependencies are managed with [pixi](https://pixi.sh). After installing pixi:

```bash
git clone https://github.com/AI-Biome/PRISTINE.git
cd PRISTINE
pixi install          # creates ./.pixi with python, minimap2, mafft, primer3-py, ...
pixi run pytest       # optional: run the test suite
```

Run the pipeline from the directory that holds your `config.yaml`:

```bash
pixi run --manifest-path /path/to/PRISTINE/pixi.toml pristine --config config.yaml
```

or activate the environment once with `pixi shell` and then call `pristine` directly.

```text
pristine --help
usage: pristine [-h] [-c CONFIG] [-t THREADS] [-o OUTPUT] [--log-level {DEBUG,INFO,WARNING,ERROR}] [--version]
```

`--threads` and `--output` override `max_cores` and `output_dir` from the config.

### Containers

A ready-made image with the locked pixi environment is on Docker Hub as
[`gaiusaugustus/pristine`](https://hub.docker.com/r/gaiusaugustus/pristine). No local
installation is needed. Run it from the directory that holds your `config.yaml`.

**Docker** (Linux, macOS, Windows). The project directory is mounted at `/data`;
`--user` keeps output files owned by you.

```bash
docker pull gaiusaugustus/pristine:0.2.0
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$(pwd):/data" gaiusaugustus/pristine:0.2.0 --config config.yaml
```

**Apptainer / Singularity** (HPC clusters). The current directory is mounted
automatically.

```bash
apptainer pull pristine.sif docker://gaiusaugustus/pristine:0.2.0
apptainer run pristine.sif --config config.yaml
```

The image is built for x86-64. On Apple Silicon Macs Docker runs it under emulation,
which is slower but works.

To build the image yourself from the repository root, use
`docker build -t pristine .` or `apptainer build pristine.sif PRISTINE/PRISTINE.def`.

## Quick start with the toy dataset

`PRISTINE/toy_dataset/` contains two *Yersinia pseudotuberculosis* genomes (targets)
and two *Yersinia pestis* genomes (non-targets) plus a ready `config.yaml`:

```bash
cd PRISTINE/toy_dataset
pixi run --manifest-path ../../pixi.toml pristine --config config.yaml

# or, without installing anything, with the container
docker run --rm --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$(pwd):/data" gaiusaugustus/pristine:0.2.0 --config config.yaml
```

With 4 threads this takes one to two minutes and writes `output/`.

## Input

### Directory layout

One folder per target group, named like the group (for a single species:
`Genus_species`). Target genome files are named `<group>_<i>.fasta`. A `non-targets/`
subfolder holds genomes of closely related species that must be distinguished from
the group. A top-level `non-targets/` folder holds non-targets used for every group.
Non-target species are identified by the file name without the trailing `_<i>`
(`Yersinia_pestis_2.fasta` is species `Yersinia_pestis`). Accepted extensions:
`.fasta .fa .fas .fna`.

```text
input/
├── Yersinia_pseudotuberculosis/
│   ├── Yersinia_pseudotuberculosis_1.fasta
│   ├── Yersinia_pseudotuberculosis_2.fasta
│   └── non-targets/
│       ├── Yersinia_enterocolitica_1.fasta
│       └── Yersinia_pestis_1.fasta
└── non-targets/
    └── Salmonella_enterica_1.fasta
```

Input files are never modified. Genomes may be draft assemblies (several contigs).

### Configuration file

All relative paths are resolved relative to the directory of the config file.

```yaml
input_dir: "data/raw_genomes/"
output_dir: "output/"
max_cores: 8                      # used for minimap2 and for parallel mafft runs
snp_avg_prop_threshold: 0.005     # loci scoring below this are not copied to informative_loci/

primer3:
  global_params:                  # any Primer3 global tag; these are the defaults
    PRIMER_NUM_RETURN: 5          # candidate pairs per locus that are checked in silico
    PRIMER_MIN_SIZE: 18
    PRIMER_OPT_SIZE: 20
    PRIMER_MAX_SIZE: 25
    PRIMER_MIN_TM: 57.0
    PRIMER_OPT_TM: 60.0
    PRIMER_MAX_TM: 63.0
    PRIMER_PRODUCT_SIZE_RANGE: "500-1800"

locus_discovery:
  window_size: 2000               # bp; must be >= the largest product size
  step: 500                       # bp between window starts on the reference genome
  minimap2_preset: asm20          # asm5 / asm10 / asm20 (increasing divergence tolerated)
  min_query_coverage: 0.9         # a window must align over >= 90 % of its length
  max_divergence: 0.15            # drop windows more divergent than this from any genome
  max_candidate_loci: 5000        # most divergent (vs non-targets) windows kept for scoring
  require_in_nontargets: true     # loci must be present in non-target genomes too
  refine_with_mafft: false        # re-align loci with mafft before primer design

snp_primer_design:
  snp_window_size: 250            # window scanned for informative-site clusters
  snp_top_n: 5                    # stop after this many loci have an accepted primer pair
  min_snps: 3                     # informative sites required per non-target species in the amplicon
  max_loci_to_try: 50
  require_snps_for_all_nontargets: true

primer_selection:
  conservation: all               # "all": primer sites identical in every genome; "targets": only in targets
  use_thermodynamics: true        # score candidate binding sites by duplex Tm (primer3)
  annealing_temp: null            # default PRIMER_MIN_TM - 3
  binding_tm_margin: 5.0          # a site binds if its duplex Tm >= annealing_temp - margin
  max_mismatches: 6          # prefilter: mismatches tolerated over the whole primer
  three_prime_region: 8      # 3'-terminal bases scrutinised separately
  three_prime_max_mismatches: 1  # mismatches tolerated inside the 3' region
  seed_length: 4             # exact k-mer used to find candidate sites
  max_offtarget_product_size: 4000
  background_dir: null            # folder of FASTA genomes in which the primers must not amplify (set this!)
amplicon_screen:
  background_identity: true       # minimap2 identity of each selected amplicon in every background genome
  remote_blast: true              # megablast the selected amplicons against NCBI nt (network, minutes; false offline)
  blast_db: nt
  hitlist_size: 100
  entrez_query: null              # e.g. "bacteria[organism]"
  min_identity_report: 90.0
make_plots: true
```

`PRIMER_MAX_NS_ACCEPTED` is always 0: the design template is N-masked at every
position that is not identical across the genomes, which is what keeps primers in
conserved sites. Config files of the previous prokka/panaroo based version are still
read; their obsolete keys are ignored with a warning.

## How it works

![PRISTINE workflow](PRISTINE-1.svg)

The figure is generated by `python docs/make_workflow_figure.py`.

1. **Genome manifest.** Genomes are discovered from the directory layout. The largest
   target genome of a group becomes the reference.
2. **Locus discovery** (minimap2 [1]). The reference is tiled into overlapping windows
   (`window_size`, `step`). All windows are mapped to every genome. A window is kept
   when, in each genome, it has exactly one alignment covering >= `min_query_coverage`
   of the window and no secondary alignment covering more than half of it. This removes
   repeats, mobile elements and multi-copy genes, which are a major source of
   unexpected PCR products. Mapping four 4.7 Mb genomes takes about ten seconds.
3. **Alignment from the mapping.** The base-level `cs` strings of the mappings are
   stacked into one reference-anchored alignment per locus. No multiple-sequence
   aligner is run here, so the number of genomes does not blow up the runtime.
   Optionally (`refine_with_mafft`), only the few loci that reach primer design are
   re-aligned with MAFFT [6].
4. **Scoring.** A column is informative for non-target species S if the bases shown by
   the target genomes and the bases shown by S's genomes are disjoint. A run of gap
   columns is one indel event. Loci are ranked by the mean proportion of informative
   sites over non-target species (and, for groups with several target species, by
   pairwise between-species sites).
5. **Primer design** (Primer3 [2, 3] through primer3-py [14]). The design template is the alignment consensus
   with every non-conserved column (and every column touched by an insertion) masked
   as N. The densest `snp_window_size` window of informative sites is passed to
   Primer3 as `SEQUENCE_TARGET`, so the amplicon must enclose it.
6. **In silico PCR.** Every candidate pair is tested against all genomes of the group.
   Binding sites are found on both strands in two stages: a NumPy [11] seed scan that finds
   every genomic window with at most `three_prime_max_mismatches` mismatches in the
   3'-terminal `three_prime_region` bases (every exact `seed_length`-mer of that region
   is looked up, so a single 3' mismatch, even at the terminal base, cannot hide a
   site) and at most `max_mismatches` overall, also at a shift of one base so template
   bulges are not missed; then Primer3's thermodynamic heterodimer calculation, based on
   nearest-neighbour parameters [4, 5], of the duplex Tm of the
   primer on that genomic window. A site binds if the Tm is at least `annealing_temp -
   binding_tm_margin` (49 C by default). All products up to
   `max_offtarget_product_size` from any primer combination are enumerated. A pair is
   accepted only if every genome gives exactly one product, located at the locus, and
   the amplicon contains at least `min_snps` informative sites for each non-target
   species.
7. **Background screen.** All pairs that survived step 6 are screened in one streaming
   pass over the genomes in `background_dir` (loaded one file at a time, so hundreds of
   genomes are fine). Any product in a background genome rejects the pair. Then the
   best surviving pair per locus is selected, one locus per genomic region, until
   `snp_top_n` loci are found. Without `background_dir` the pipeline warns loudly:
   the screen then knows nothing beyond your input genomes.
8. **Amplicon screen.** Primer screening asks *what amplifies*; this step asks *who else
   carries this marker*. For every selected locus the reference amplicon is (a) mapped
   to each background genome with minimap2 and (b), unless `remote_blast: false`,
   megablasted [7, 8] against NCBI nt [9] through NCBI's web service, submitted and
   parsed with Biopython [10]. Each hit is labelled target,
   known non-target or unknown taxon by its organism name. An unknown organism whose
   identity to the amplicon is at least the *tolerated identity* (100 % minus the
   fraction of informative sites separating the target from its nearest designed
   non-target) is a false-positive risk and is reported with a warning: add it as a
   non-target or background genome and rerun.

## Output

```text
output/
├── run_summary.csv                      one row per group
├── pristine.log
└── <group>/
    ├── locus_discovery/                 reference windows and minimap2 PAF files
    ├── loci/
    │   ├── locus_table.csv              locus -> coordinates in every genome
    │   └── aligned/<locus>.fasta        reference-anchored alignments of all candidate loci
    ├── informative_loci/
    │   ├── snp_summary.csv              scores and informative positions per locus (see below)
    │   ├── heatmap_informativeness.png
    │   ├── <locus>.fasta                alignments passing snp_avg_prop_threshold
    │   ├── consensus_sequences/         IUPAC consensus and N-masked design template per designed locus
    │   └── snp_density_plots/           informative-site histograms of selected loci
    └── primers/
        ├── selected_primers.csv         the accepted primer pair per locus  <-- main result
        ├── primer_design_summary.csv    every candidate pair with pass/fail reasons
        ├── locus_design_log.csv         per locus: target region and outcome
        ├── <locus>_primers.fasta        LEFT/RIGHT sequences of accepted pairs
        ├── amplicons/<locus>_amplicons.fasta  predicted amplicon of every genome, locus orientation
        └── amplicon_screen/
            ├── amplicon_screen_summary.csv  per locus: tolerated identity, closest background genome, risk taxa
            ├── background_identity.csv      amplicon vs. background genome identities
            ├── blast_hits.csv               NCBI hits with organism, identity, status (remote_blast only)
            └── reference_amplicons.fasta    query sequences, ready for manual Primer-BLAST / BLAST
```

Locus names are `<contig>_<start>_<end>` on the reference genome.

`snp_summary.csv` columns: `Combined_Informativeness_Score` (ranking),
`Avg_Prop_Informative_SNPs_External`/`_Internal`, and per non-target species `S`:
`Prop_External_S`, `Abs_External_S` (substitutions + indel events),
`Subst_External_S`, `Indel_External_S`, `SNP_Pos_External_S` (alignment columns) and
`Core_Region_External_S`. Groups with several target species get the same columns as
`..._Internal_<A>_vs_<B>`.

`selected_primers.csv` / `primer_design_summary.csv` columns include the primer
sequences, Tm and GC, `Template_Product_Size` (on the alignment, including insertion
columns of other genomes), `On_Target_Product_Sizes` (real predicted size per genome),
`Sites_In_Amplicon_External_S`, `Min_On_Target_Site_Tm` (weakest predicted primer
duplex among the intended sites) and `Failure_Reasons` for rejected pairs.

## Notes on specificity

Primers with Tm around 60 C, no ambiguous bases, and a single predicted product in
every input genome are necessary but not sufficient for a clean PCR: the in silico
check only knows the genomes you provide. Add genomes of organisms expected in your
samples (host, common contaminants, other community members, e.g. the taxa seen in
16S profiles of the same samples) to `background_dir`, or as non-targets if they
should also be discriminated by the amplicon. For a final shortlist, Primer-BLAST
against NCBI nt remains the broadest available screen. The thermodynamic binding model
is a prediction: a reported duplex Tm just above the threshold means a weak but
possible product, not a certain one.

## Development

```bash
pixi run test                  # pytest on synthetic genomes (uses minimap2, mafft, primer3)
```

The code lives in `src/pristine/`: `manifest` (input discovery), `loci` (windows,
minimap2), `pseudoalign` (alignments from cs tags), `scoring`, `consensus`, `primers`
(primer3), `insilico_pcr`, `pipeline`, `cli`.

## References

If you use PRISTINE, please also cite the tools it builds on.

1. Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics* 34(18):3094-3100. <https://doi.org/10.1093/bioinformatics/bty191>
2. Untergasser A. et al. (2012). Primer3: new capabilities and interfaces. *Nucleic Acids Research* 40(15):e115. <https://doi.org/10.1093/nar/gks596>
3. Koressaar T., Remm M. (2007). Enhancements and modifications of primer design program Primer3. *Bioinformatics* 23(10):1289-1291. <https://doi.org/10.1093/bioinformatics/btm091>
4. SantaLucia J. Jr. (1998). A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95(4):1460-1465. <https://doi.org/10.1073/pnas.95.4.1460>
5. SantaLucia J. Jr., Hicks D. (2004). The thermodynamics of DNA structural motifs. *Annual Review of Biophysics and Biomolecular Structure* 33:415-440. <https://doi.org/10.1146/annurev.biophys.32.110601.141800>
6. Katoh K., Standley D.M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. *Molecular Biology and Evolution* 30(4):772-780. <https://doi.org/10.1093/molbev/mst010>
7. Camacho C. et al. (2009). BLAST+: architecture and applications. *BMC Bioinformatics* 10:421. <https://doi.org/10.1186/1471-2105-10-421>
8. Zhang Z., Schwartz S., Wagner L., Miller W. (2000). A greedy algorithm for aligning DNA sequences. *Journal of Computational Biology* 7(1-2):203-214. <https://doi.org/10.1089/10665270050081478>
9. Sayers E.W. et al. (2026). Database resources of the National Center for Biotechnology Information in 2026. *Nucleic Acids Research* 54(D1):D20-D27. <https://doi.org/10.1093/nar/gkaf1060>
10. Cock P.J.A. et al. (2009). Biopython: freely available Python tools for computational molecular biology and bioinformatics. *Bioinformatics* 25(11):1422-1423. <https://doi.org/10.1093/bioinformatics/btp163>
11. Harris C.R. et al. (2020). Array programming with NumPy. *Nature* 585(7825):357-362. <https://doi.org/10.1038/s41586-020-2649-2>
12. McKinney W. (2010). Data Structures for Statistical Computing in Python. *Proceedings of the 9th Python in Science Conference (SciPy 2010)*. <https://doi.org/10.25080/Majora-92bf1922-00a>
13. Hunter J.D. (2007). Matplotlib: A 2D Graphics Environment. *Computing in Science & Engineering* 9(3):90-95. <https://doi.org/10.1109/MCSE.2007.55>
14. primer3-py: Python bindings for Primer3. <https://github.com/libnano/primer3-py>

Tables are written with pandas [12] and plots drawn with matplotlib [13]. The environment
is managed with [pixi](https://pixi.sh); containers use [Docker](https://www.docker.com)
or [Apptainer](https://apptainer.org).
