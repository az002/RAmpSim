# RAmpSim - Read Amplification Simulator

A thermodynamically-aware simulator for bait-capture sequencing experiments. RAmpSim generates realistic DNA fragments based on probe alignments, thermodynamic binding calculations, and genomic abundance profiles.

## Overview

RAmpSim simulates long-read sequencing data for experiments where DNA fragments are generated through a combination of:
- **Bait-capture** (e.g., probes targeting antibiotic resistance genes or other regions of interest)
- **Random background amplification**

The simulator uses nearest-neighbor thermodynamic parameters to calculate probe-target binding affinities and generates fragment length distributions following a log-normal model, mimicking real sequencing data characteristics.

## Installation

### Prerequisites

- Rust toolchain (stable or nightly)
- Cargo package manager

### Build from Source

```bash
git clone https://github.com/az002/RAmpSim.git
cd RAmpSim
cargo build --release
```

The compiled binary will be available at `target/release/read_sim`.

## Quick Start

```bash
./target/release/read_sim \
  --ref-path references.fa \
  --sam-path probe_alignments.sam \
  --out-path simulated_reads.fa \
  --hits-path hits.txt \
  --abundances abundances.txt \
  --seqids sequence_mapping.txt \
  --temperature 60.0 \
  --nfrag 10000 \
  --split 0.7
```

This generates 10,000 fragments: 70% from probe-enriched targets, 30% from random background.

## Usage

### Command Line Arguments

#### Required Arguments

| Argument | Short | Description |
|----------|-------|-------------|
| `--ref-path` | `-r` | Reference sequences (FASTA format) |
| `--sam-path` | `-s` | Probe alignments (SAM/BAM format) |
| `--out-path` | `-o` | Output file for simulated reads (FASTA) |
| `--abundances` | `-a` | Tab-separated file with reference abundances |
| `--seqids` | | Mapping of sequence IDs to reference IDs |
| `--temperature` | `-t` | Hybridization temperature (°C) |
| `--nfrag` | `-n` | Total number of fragments to generate |
| `--split` | | Fraction of probe-enriched fragments (0.0-1.0) |

#### Optional Arguments

| Argument | Short | Default | Description |
|----------|-------|---------|-------------|
| `--flen` | `-l` | 8000 | Mean fragment length (bp) |
| `--lognorm-sd` | | 0.4 | Standard deviation for log-normal distribution |
| `--mult-file` | | None | Probe multiplicity file (tab-separated) |
| `--hits-path` | | None | Output path for probe binding site scores |

### Input File Formats

#### 1. Reference Sequences (`--ref-path`)

Standard FASTA format:
```fasta
>sequence_id_1
ATCGATCGATCGATCG...
>sequence_id_2
GCTAGCTAGCTAGCTA...
```

#### 2. Probe Alignments (`--sam-path`)

SAM or BAM format with probe sequences aligned to references. Can be generated with tools like:
- Bowtie2
- BWA-MEM
- Minimap2

```sam
@HD	VN:1.0	SO:unsorted
@SQ	SN:chr1	LN:248956422
probe_001	0	chr1	1000	60	50M	*	0	0	ATCGATCG...	IIIIIIII...
```

#### 3. Abundances File (`--abundances`)

Tab-separated: reference_id → relative abundance
```
E_coli_K12	0.45
S_aureus	0.30
P_aeruginosa	0.25
```

**Note**: Abundances don't need to sum to 1.0; they represent relative proportions.

#### 4. Sequence ID Mapping (`--seqids`)

Tab-separated: sequence_id → reference_id
```
NC_000913.3	E_coli_K12
NC_007795.1	S_aureus
NC_002516.2	P_aeruginosa
```

Maps individual sequences (chromosomes, plasmids, contigs) to reference genomes.

#### 5. Probe Multiplicity File (Optional, `--mult-file`)

Tab-separated: probe_name → copy_number
```
telomere_probe_1	10
telomere_probe_2	5
repeat_element_1	3
```

Higher values increase sampling probability for that probe.

### Output Files

#### 1. Simulated Reads FASTA (`--out-path`)

**Probe-enriched fragments:**
```fasta
>probe_name|sequence_id|probe_start|fragment_start-fragment_end
ATCGATCGATCGATCG...
```

**Background fragments:**
```fasta
>Background|sequence_id|fragment_start-fragment_end
GCTAGCTAGCTAGCTA...
```

#### 2. Hits File (Optional, `--hits-path`)

Tab-separated file with probe alignment scores. Only generated if `--hits-path` is specified.

```
probe_001	NC_000913.3	1000	1050	0.234
probe_001	NC_007795.1	5000	5050	0.156
probe_002	NC_000913.3	2000	2050	0.423
```

Columns:
1. Probe name
2. Sequence ID where probe aligned
3. Alignment start position
4. Alignment end position
5. Normalized probability score (based on thermodynamics × abundance)

#### 3. Standard Output Logging

Real-time fragment generation log:
```
[TELSEQ]:E_coli_K12:NC_000913.3:15234-23451:probe_1
[TELSEQ]:S_aureus:NC_007795.1:8923-17234:probe_2
[BG]:E_coli_K12:NC_000913.3:45678-53890:
[BG]:P_aeruginosa:NC_002516.2:12345-20567:
```

Format: `[TYPE]:reference_id:sequence_id:start-end:probe_name`

## Examples

Example files are provided in the examples directory. Download reference genomes with
```
wget https://s3.amazonaws.com/zymo-files/BioPool/ZymoBIOMICS.STD.refseq.v2.zip
unzip ZymoBIOMICS.STD.refseq.v2.zip
cat ZymoBIOMICS.STD.refseq.v2/Genomes/*.fasta > references.fa
```
