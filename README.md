> [!WARNING]
> This branch is based on the old prokka/panaroo version of PRISTINE. `main` was
> rewritten for PRISTINE 0.2.0 (package in `src/pristine/`, minimap2 locus discovery,
> thermodynamic in silico PCR, pixi environment), so this branch can no longer be
> merged into `main` without porting the changes by hand. Use `main` for new work.

# Marker Loci Identification and Primer Design Tool

This tool aims to find a minimal set of reliable markers for identification of given target species. Optionally, for each target species, a list of non-target closely related species may be specified.

The marker regions are searched for within single-copy core genes only.

1. Pan-genome analysis

Input genomes are annotated using Prokka (Seemann, 2014) and core genes are identified with Panaroo (Tonkin-Hill, 2020). Annotated sequences are saved to a local PostgreSQL database. Single-copy core genes are then queried from this database and saved in the multi-FASTA format.

2. Multiple sequence alignment

Single copy core genes are aligned using one of the following algorithms:

* ClustalW
* MUSCLE
* MAFFT
* probabilistic multiple alignment program PRANK (Löytynoja, 2014)

3. Identification of conserved regions

Identification of conserved regions for each species can be carried out by means of:
* Consensus sequence
* Shannon's entropy
* Quasi-alignment-based method proposed by Nagar & Hahsler (2013)

After finding species-specific marker regions, a coverage matrix is created showing which of the found markers can be used to reliably identify which of the target isolates. Based on this matrix, a minimal set of marker regions is selected.

## Requirements

- Python 3.x
- External tools: Prokka, Panaroo, ClustalW, MUSCLE, MAFFT, PRANK, Minimap2, Primer3 (if using default primer design algorithm)
- A reference database compatible with Minimap2 for sequence alignment

## Usage

```bash
python marker_primer_tool.py [options]

# Options

Marker Loci Identification
-a, --alignment-tool (required): Alignment tool to use. Possible values: clustalw, muscle, mafft, prank.
-r, --reference-db (required): Path to the reference database for Minimap2.
-I, --identity-threshold: Identity threshold for Minimap2 (0-100).
-P, --proportion-threshold: Threshold for the proportion of correctly identified sequences (0-100).
-l, --min-length: Minimum length of conserved regions.
-g, --algorithm: Algorithm for finding conserved regions. Possible values: consensus_sequence, shannon_entropy, quasi_alignment.

Primer Design
-i, --input_file (required): The input file to process.
-o, --output_file (default: output.txt): The output file to write to.
-s, --strategy (required): Choose a general strategy - primer design only for the WGS dataset (amplicon) or identification of marker loci in the WGS dataset followed by the primer design for the identified marker loci (marker).
-a, --primer_design_algorithm (required): Choose a primer design algorithm. Possible values: primer3, custom.
-p, --primer3_parameters: Path to the Primer3 parameters config file.
-S, --primer_summarizing_algorithm (required): An algorithm for summarizing primers. Possible values: frequency, consensus.
-c, --specificity_check_algorithm (required): An algorithm for primer specificity checking. Only blast is currently supported.
-d, --database_for_specificity_check (required): A database for checking the specificity of primers.
-n, --n_most_frequent: A number of the most frequently occurring primers to further work with.