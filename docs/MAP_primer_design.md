# Designing MAP primers with PRISTINE

How to design ONT amplicon primers that identify *Mycobacterium avium* subsp.
*paratuberculosis* (MAP) and tell it apart from the other *M. avium* subspecies.
It covers the genome data to use, how to download it, the parameters to set and
how to run the pipeline and read the results.

## 1. Design decisions

| Role in PRISTINE | Taxa | Why |
|---|---|---|
| Target group | MAP, C-type (cattle) and S-type (sheep) strains | The primers must amplify every MAP lineage. A primer site missing in S-type strains means false negatives in sheep and goat samples. |
| Non-targets (amplified, and told apart by sequence) | *M. avium* subsp. *hominissuis* (MAH), subsp. *avium* (MAA), subsp. *silvaticum* (MAS) | These are MAP's closest relatives and are common in feces, soil and water. MAP vs MAA is about 98.5 % ANI (skani, MAP E1 vs MAA DSM 44156), so 1800 bp amplicons can carry enough informative sites. |
| Background (must not amplify) | *M. intracellulare*, *M. chimaera*, other MAC species, other mycobacteria (including *M. bovis*), high-GC actinobacteria, host genome | Any predicted product in these genomes rejects a primer pair. |

The classical MAP PCR targets (IS900, F57, ISMav2) are MAP-specific insertions and
mostly multi-copy. PRISTINE finds something different: single-copy loci that are
present in all subspecies and carry MAP-specific SNPs and indels. Windows that hit
several places in a genome, such as IS900 copies, are removed on purpose. You should
not expect IS900-like markers from this run, and none of its results replace them.

Put *M. intracellulare* and *M. chimaera* in the background first (section 3). If
almost every candidate pair is then rejected by the background screen, the conserved primer sites also bind in *M. intracellulare*. You
can see this in `locus_design_log.csv` as `all passing pair(s) rejected by the
background screen`, and in the `background Mycobacterium_intracellulare_*` entries of
`Failure_Reasons`. In that case, move both species into `non-targets/` and rerun. The
amplicon is then allowed to amplify them but must carry informative sites against them.

## 2. Requirements

- A PRISTINE checkout with `pixi install` done, and `pixi` on your `PATH`.
- Network access to NCBI.
- About 5 GB of free disk space while the cattle genome is downloaded and split, and
  about 4 GB after cleanup (genomes plus output).
- The NCBI `datasets` CLI and `skani` do not need to be installed. They run through
  `pixi exec`.

## 3. Download the genomes

Everything below runs in its own work directory, outside the repository.

```bash
export WORK=$HOME/pristine_map
mkdir -p $WORK/{lists,dl,qc} && cd $WORK

ds() { pixi exec -s ncbi-datasets-cli "$@"; }

# complete RefSeq genomes of a taxon -> accession, organism, strain, length
lst() {
  ds datasets summary genome taxon "$1" --assembly-source RefSeq --assembly-level complete \
     --exclude-atypical --as-json-lines |
  ds dataformat tsv genome --elide-header \
     --fields accession,organism-name,organism-infraspecific-strain,assmstats-total-sequence-len
}

# fetch <list.tsv> <dest_dir> <Label>  ->  <dest_dir>/<Label>_<i>.fasta  (+ row in genome_map.tsv)
fetch() {
  local list=$1 dest=$2 label=$3 i
  i=$(ls "$dest"/"${label}"_*.fasta 2>/dev/null | wc -l)
  mkdir -p "$dest" "dl/$label"
  cut -f1 "$list" > "dl/$label/acc.txt"
  ds datasets download genome accession --inputfile "dl/$label/acc.txt" --include genome \
     --filename "dl/$label.zip" --no-progressbar
  pixi exec -s unzip unzip -q -o "dl/$label.zip" -d "dl/$label"
  for f in dl/"$label"/ncbi_dataset/data/GC*/*.fna; do
    i=$((i+1))
    cp "$f" "$dest/${label}_$i.fasta"
    printf '%s\t%s\n' "$(basename "$(dirname "$f")")" "${label}_$i" >> genome_map.tsv
  done
  rm -rf "dl/$label" "dl/$label.zip"
}
```

### 3.1 List the candidates and review them

```bash
lst "Mycobacterium avium subsp. paratuberculosis" > lists/MAP.tsv
lst "Mycobacterium avium subsp. hominissuis"      > lists/MAH.tsv
lst "Mycobacterium avium subsp. avium"            > lists/MAA.tsv
lst "Mycobacterium avium subsp. silvaticum"       > lists/MAS.tsv
lst "Mycobacterium intracellulare"                > lists/MINT_all.tsv   # includes subsp. chimaera
grep -v chimaera lists/MINT_all.tsv > lists/MINT.tsv
grep    chimaera lists/MINT_all.tsv > lists/MCHIM.tsv
wc -l lists/*.tsv
```

On 2026-11-07, NCBI RefSeq had 28 complete MAP genomes, 42 MAH, 4 MAA, 1 MAS and 42
*M. intracellulare* (11 of them subsp. *chimaera*).

Review `lists/MAP.tsv` before downloading:

- Keep both MAP lineages. Check the strain names for S-type (sheep) isolates. If
  the complete genomes contain none, add S-type genomes from a lower assembly level
  (repeat `lst` with `--assembly-level complete,chromosome,scaffold` and pick them by
  strain name). Draft assemblies are accepted. Their only cost is that loci split
  across contigs are dropped.
- K-10 is the MAP type genome and appears
  twice (`GCF_000007865.1` and the PacBio re-sequencing `GCF_025665475.1`). Keeping
  both does no harm.
- Don't remove a target genome to rescue a design. Each one you drop is a MAP
  lineage the primers are no longer tested on.

### 3.2 Targets and non-targets

The folder name sets the group name. File names give the species label: PRISTINE
strips the trailing `_<i>`, so `Mycobacterium_avium_hominissuis_7.fasta` is counted as
species `Mycobacterium_avium_hominissuis`.

```bash
G=input/Mycobacterium_avium_paratuberculosis
fetch lists/MAP.tsv $G             Mycobacterium_avium_paratuberculosis
fetch lists/MAH.tsv $G/non-targets Mycobacterium_avium_hominissuis
fetch lists/MAA.tsv $G/non-targets Mycobacterium_avium_avium
fetch lists/MAS.tsv $G/non-targets Mycobacterium_avium_silvaticum
```

Do not add genomes labelled only "*Mycobacterium avium*", with no subspecies. Most of
them are MAH, but the label is not reliable. Add one only after the ANI check in 3.4
places it, and then name it after its subspecies.

### 3.3 Background genomes

`background_dir` is read non-recursively: every FASTA file directly inside
`background/` is one background genome.

```bash
fetch lists/MINT.tsv  background Mycobacterium_intracellulare
fetch lists/MCHIM.tsv background Mycobacterium_chimaera

cat > lists/background_other.tsv <<'EOF'
GCF_000222105.3	Mycobacterium_colombiense
GCF_002285715.1	Mycobacterium_marseillense
GCF_000195955.2	Mycobacterium_tuberculosis
GCF_000195835.3	Mycobacterium_bovis
GCF_000157895.3	Mycobacterium_kansasii
GCF_002101675.1	Mycobacterium_gordonae
GCF_002086735.1	Mycobacterium_scrofulaceum
GCF_000069185.1	Mycobacteroides_abscessus
GCF_000015005.1	Mycolicibacterium_smegmatis
GCF_003013675.1	Rhodococcus_equi
GCF_000009805.1	Nocardia_farcinica
GCF_000203835.1	Streptomyces_coelicolor
GCF_003932295.1	Corynebacterium_bovis
GCF_006739385.1	Cutibacterium_acnes
GCF_002706665.1	Bifidobacterium_pseudolongum
GCF_000005845.2	Escherichia_coli
EOF
while read -r acc label; do
  printf '%s\n' "$acc" > lists/one.tsv
  fetch lists/one.tsv background "$label"
done < lists/background_other.tsv
```

All accessions above were checked against NCBI on 2026-11-07. The list is a
starting point. Add any taxa you find in 16S profiles of the actual samples (fecal
or milk community, kit and lab contaminants).

The host genome needs splitting. The in silico PCR loads one background file at a
time and builds a k-mer index of it in memory. A whole 2.8 Gb mammalian genome in one file would need
tens of GB of RAM, so split it into one file per chromosome. The largest bovine
chromosome then needs about 3 GB (measured on a 160 Mb test chromosome).

Background genomes are screened in file-name order, and a pair rejected in one genome
is not screened in the rest. The host files are therefore named `host_*`: lower-case
names sort after the capitalised species files, so the close relatives (*M.
intracellulare* etc.), which reject most pairs, are screened first and the slow
mammalian chromosomes only see the few pairs that are left. Watch the
`background i/N screened` lines in `output/pristine.log`.

```bash
ds datasets download genome taxon "Bos taurus" --reference --include genome \
   --filename dl/bos.zip --no-progressbar             # GCF_002263795.3, ~2.8 Gb
pixi exec -s unzip unzip -q -o dl/bos.zip -d dl/bos
awk '/^>/{ id=substr($1,2)
           f = (id ~ /^NC_/) ? "background/host_Bos_taurus_" id ".fasta" : "background/host_Bos_taurus_unplaced.fasta" }
     { print > f }' dl/bos/ncbi_dataset/data/GCF_*/*.fna
rm -rf dl/bos dl/bos.zip
```

For sheep or goat samples, do the same with `"Ovis aries"` (GCF_016772045.2) or
`"Capra hircus"` (GCF_001704415.2, which has about 30,000 unplaced scaffolds that all
go into the `unplaced` file).

### 3.4 Check species labels with ANI

NCBI subspecies labels in MAC are sometimes wrong. A mislabelled genome breaks the
design in one of two ways. A MAH genome filed as MAP removes every informative site
against MAH. A MAP genome filed as MAH makes the target and non-target sets overlap.

```bash
pixi exec -c conda-forge -c bioconda -s skani \
  skani triangle -t 6 -E $G/*.fasta $G/non-targets/*.fasta > qc/ani.tsv 2>/dev/null

# min / max ANI for each pair of species labels
awk -F'\t' 'NR>1{ n=split($1,a,"/"); m=split($2,b,"/"); sa=a[n]; sb=b[m]
     sub(/_[0-9]+\.fasta$/,"",sa); sub(/_[0-9]+\.fasta$/,"",sb)
     k=(sa<sb)? sa" vs "sb : sb" vs "sa
     if(!(k in mn)||$3<mn[k]) mn[k]=$3; if(!(k in mx)||$3>mx[k]) mx[k]=$3 }
     END{ for(k in mn) printf "%-75s min %6.2f  max %6.2f\n", k, mn[k], mx[k] }' qc/ani.tsv | sort
```

What to look for:

- The MAP vs MAP minimum should be clearly higher than the MAP vs MAH and MAP vs MAA
  maximums. MAP is near-clonal.
- If the maximum of a MAP-vs-other pair reaches the MAP-vs-MAP range, that genome is
  probably mislabelled. Find it with
  `awk -F'\t' '$3>99.5' qc/ani.tsv | grep paratuberculosis | grep non-targets`.
  Move or delete the file and look up its accession in `genome_map.tsv`.

## 4. Configuration

Save this as `$WORK/config.yaml`. Relative paths are resolved from the config file's
directory. Comments mark where it differs from the PRISTINE defaults and why.

```yaml
input_dir: "input/"
output_dir: "output/"
max_cores: 6                      # leave a core or two free for the rest of the system
snp_avg_prop_threshold: 0.005

primer3:
  global_params:
    PRIMER_NUM_RETURN: 5          # candidate pairs per locus; each one costs an in silico PCR over all genomes
    PRIMER_MIN_SIZE: 18
    PRIMER_OPT_SIZE: 20
    PRIMER_MAX_SIZE: 25
    PRIMER_MIN_TM: 58.0           # MAP is ~69 % GC: a slightly higher Tm window
    PRIMER_OPT_TM: 62.0           #   finds more acceptable conserved sites
    PRIMER_MAX_TM: 65.0
    PRIMER_PAIR_MAX_DIFF_TM: 3.0
    PRIMER_MIN_GC: 45.0
    PRIMER_MAX_GC: 70.0           # default 65 rejects many primers in a 69 % GC genome
    PRIMER_MAX_POLY_X: 4          # keeps G/C runs short
    PRIMER_PRODUCT_SIZE_RANGE: "800-1800"   # longer amplicons carry more informative sites

locus_discovery:
  window_size: 2000               # >= largest product size
  step: 500
  minimap2_preset: asm20
  min_query_coverage: 0.9
  max_divergence: 0.15
  max_candidate_loci: 2000        # default 5000; only the top loci are tried, and each
                                  # candidate writes a ~90-genome alignment to disk
  require_in_nontargets: true
  refine_with_mafft: false

snp_primer_design:
  snp_window_size: 250
  snp_top_n: 5
  min_snps: 3                     # minimum sites per non-target subspecies; more is better for ONT reads
  max_loci_to_try: 100            # default 50; MAA/MAS are close to MAP, so more loci fail
  require_snps_for_all_nontargets: true

primer_selection:
  conservation: all
  use_thermodynamics: true
  annealing_temp: null            # = PRIMER_MIN_TM - 3 = 55 C
  binding_tm_margin: 5.0          # off-target sites with duplex Tm >= 50 C count as binding
  max_mismatches: 6
  three_prime_region: 8
  three_prime_max_mismatches: 1
  seed_length: 4
  max_offtarget_product_size: 4000
  background_dir: "background/"

amplicon_screen:
  background_identity: true
  remote_blast: true              # megablast of the selected amplicons against NCBI nt (a few minutes)
  blast_db: nt
  hitlist_size: 100
  entrez_query: null              # whole nt, so hits outside bacteria show up too
  min_identity_report: 90.0

make_plots: true
```

If a run gives fewer than `snp_top_n` loci, relax in this order and rerun:

1. `max_loci_to_try: 200` and `max_candidate_loci: 4000`.
2. `PRIMER_PRODUCT_SIZE_RANGE: "500-1800"`.
3. `min_snps: 2`, the last resort. With fewer sites, single ONT reads become harder to
   assign to a subspecies.

Do not set `conservation: targets` or `require_snps_for_all_nontargets: false` for this
use case. With `conservation: targets`, MAH could fail to amplify, so a negative
result could no longer tell "MAH present" from "nothing there". With
`require_snps_for_all_nontargets: false`, a marker could be blind to MAA.

## 5. Run

```bash
cd $WORK
export PRISTINE_DIR=/path/to/PRISTINE   # your checkout
P=$PRISTINE_DIR/pixi.toml

# 1. check the installation (one to two minutes)
(cd $PRISTINE_DIR/PRISTINE/toy_dataset && pixi run --manifest-path $P pristine --config config.yaml --threads 4)

# 2. design run
pixi run --manifest-path $P pristine --config config.yaml --threads 6 2>&1 | tee run.log
```

To get a quick yield check before the full screen, run once with
`background_dir: null` and `remote_blast: false`. Then rerun with both turned on.

Most of the runtime goes into the in silico PCR of each candidate pair against every
input genome. It runs in `max_cores` worker processes, one pair per worker, on Linux;
on macOS it runs on one core. Measured on 14 real genomes (8 MAP, 6 MAH; 8-thread
machine, `max_loci_to_try: 30`, no background): 268 s on 1 core, 110 s with
`--threads 6`, identical results. With the full set (~75 input genomes) one pair takes
about 30 s of CPU time, so a design that tests a few hundred pairs needs a few CPU
hours, divided by the number of cores. More cores help directly here (16 threads
roughly halve the time of 8). The background screen adds the host chromosomes (about
20 min to index the cattle genome, then a few minutes per surviving pair, shared over
the cores) and the remote BLAST.

## 6. Read the results

All paths are below `output/Mycobacterium_avium_paratuberculosis/`.

| File | What to check |
|---|---|
| `primers/selected_primers.csv` | The main result. For each pair, `Sites_In_Amplicon_External_<subsp>` should be at least 3 for MAH, MAA and MAS. Higher numbers make reads easier to classify. `On_Target_Product_Sizes` should be similar across all MAP genomes. A low `Min_On_Target_Site_Tm` means a weak primer site in some genome. |
| `primers/locus_design_log.csv` | Why loci failed. If many say `rejected by the background screen`, apply the fallback in section 1. |
| `primers/primer_design_summary.csv` | Every candidate pair with its `Failure_Reasons`. |
| `primers/amplicon_screen/amplicon_screen_summary.csv` | Organisms with a sequence close to the amplicon that are not in your input. Any `BLAST_Risk_Taxa` entry has to be explained or added as a non-target or background genome, followed by a rerun. |
| `primers/amplicon_screen/background_identity.csv` | How similar each amplicon is to the background genomes. A background genome at or above the tolerated identity could be mistaken for MAP from the reads. |
| `informative_loci/snp_density_plots/` | Where the informative sites sit in each amplicon. Evenly spread is better than one cluster. |

Before ordering primers:

1. Run Primer-BLAST on the shortlisted pairs against nt, and against the *Bos taurus*
   RefSeq genome if the host was not in `background/`. The sequences are in
   `primers/<locus>_primers.fasta`.
2. Order the top 2–3 pairs, each from a different locus.
3. Wet-lab panel: MAP C-type and S-type DNA, MAH, MAA, *M. intracellulare*, a MAP-free
   bovine fecal DNA extract, and a no-template control. Run a gradient PCR around the
   55 °C design annealing temperature. The template is high-GC, so test DMSO or
   betaine if the yield is low.
