# 🧬 GBS Pipeline: Comprehensive QC Pipeline for *Streptococcus agalactiae* WGS

**GBS Pipeline** is a Nextflow DSL2 pipeline that performs comprehensive quality control on *Streptococcus agalactiae* (Group B Streptococcus, GBS) whole-genome sequencing (WGS) data.  
It includes read-level, taxonomic, assembly, mapping, and variant QC, and generates structured reports for downstream use.

> 🔧 Based on the [Global PneumoSeq (GPS) pipeline](https://github.com/GlobalPneumoSeq/gps-pipeline), GBS Pipeline is adapted to focus specifically on high-quality, assembly-backed, taxonomically verified *S. agalactiae* genomes.

---

## 🔍 What the Pipeline Does

The pipeline executes the following QC stages **for every sample**:

| QC Stage         | Tool(s) Used                            | Purpose                                        |
|------------------|------------------------------------------|------------------------------------------------|
| **FASTQ check**  | Internal validation script               | Ensures correct naming and paired-end layout   |
| **Read QC**      | `fastp`                                  | Adapter trimming and per-read quality check    |
| **Assembly**     | `shovill` (default) or `unicycler`       | De novo genome assembly                        |
| **Assembly QC**  | `quast`                                  | Checks genome size and contig count          |
| **Mapping QC**   | `bwa`, `samtools`                        | Align reads to reference, compute depth/coverage|
| **Taxonomy QC**  | `kraken2`, `bracken`                     | Confirm species identity and screen contaminants|
| **Final Report** | Internal summary script                  | Aggregates QC metrics with pass/fail flags     |

---

## 🚀 How to Run

### 🛠️ 1. Clone the repository

```bash
git clone https://github.com/yourusername/gbs-pipeline.git
cd gbs-pipeline
```

### 📁 2. Prepare input data

1. Create an `input/` directory:
   ```bash
   mkdir input
   ```

2. Place your FASTQ files inside. Supported formats:
   - Paired-end: `SAMPLE_R1.fastq.gz` / `SAMPLE_R2.fastq.gz`
   - Accepts `_R1`, `_1`, or `_R1_001` patterns

> ⚠️ If `input/` is missing or empty, you’ll get:
```
ERROR ~ No files match pattern `*_{,R}{1,2}{,_001}.{fq,fastq}{,.gz}` at path: ./input/
```

---

## ▶️ Run Commands

### Docker-based run (recommended)

```bash
nextflow run main.nf -profile standard
```

### Singularity-based run

```bash
nextflow run main.nf -profile singularity
```

### Run on LSF cluster

```bash
nextflow run main.nf -profile lsf
```

---

## ⚙️ Key Parameters (with defaults)

| Parameter             | Description                                     | Default                          |
|------------------------|-------------------------------------------------|----------------------------------|
| `--input`             | Input FASTQ path                                | `input/*.fastq.gz`               |
| `--output`            | Output directory                                | `output/`                        |
| `--ref_genome`        | Reference genome FASTA (for BWA mapping)        | `data/CP129876.1.fasta`          |
| `--assembler`         | Assembler: `shovill` or `unicycler`             | `shovill`                        |
| `--min_contig_length` | Minimum contig length to keep in assemblies     | `500`                            |
| `--read_len`          | Read length for Bracken                         | `150`                            |


---

## ✅ Default QC Thresholds

| QC Metric                          | Threshold         |
|-----------------------------------|-------------------|
| *S. agalactiae* abundance         | ≥ 70.00%          |
| Top non-agalactiae species        | ≤ 5.00%           |
| Het-SNPs (<50-bp distance)        | < 40              |
| Reference genome coverage         | ≥ 70.0%           |
| Mean read depth                   | ≥ 20×             |
| Total contigs                     | ≤ 500             |
| Assembly genome length            | 1.4 – 2.8 Mb       |

Thresholds are set in `nextflow.config` and can be adjusted under `params`.

---

## 📦 Output Structure

Each run produces:

```
output/
├── assemblies
│   ├── sample.contigs.fasta
│  
├── sample_reports
    ├── sample_report.csv
```

---

## 🧱 Requirements

- [Nextflow](https://www.nextflow.io/) ≥ 22.04
- Docker or Singularity (for containerised runs)
- Tools (automatically pulled via containers):
  - fastp, shovill, unicycler, bwa, samtools, bcftools, kraken2, bracken, quast

---

## 📜 Licence

This pipeline is distributed under the [GNU General Public License v3.0](LICENSE),
in accordance with the licence of the original [Global PneumoSeq (GPS) pipeline](https://github.com/GlobalPneumoSeq/gps-pipeline)

---

## 🙏 Acknowledgements

This pipeline is adapted from the [Global PneumoSeq (GPS) pipeline](https://github.com/GlobalPneumoSeq/gps-pipeline), with modifications for high-quality GBS-specific WGS QC workflows.
