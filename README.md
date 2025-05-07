# microDNA Detector

This project implements a Python-based algorithm to detect **microDNA** (extrachromosomal circular DNA) in a cancer cell line from aligned sequencing data. The algorithm identifies and scores candidate circular DNA junctions using soft-clipped reads from a BAM file and validates them using microhomology detection.

---

## Project Structure

```
├── bam_read.py                      # Main microDNA detection script
├── detected_microdna_circles.txt   # Output: list of predicted circles
└── Data/
    ├── SRR413984.sorted.NC_000001.10.bam          # Input BAM file
    ├── SRR413984.sorted.NC_000001.10.bam.bai      # BAM index file
    └── GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna  # Chromosome 1 FASTA
```

---

## Requirements

This script is intended to run in a **Python 3 WSL (Linux)** environment.

Install the required packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pysam parasail tqdm numpy
sudo apt install samtools  # for indexing FASTA
```

---

## Setup

Ensure the following:

1. Your BAM file is sorted and indexed (`.bam` + `.bai`).
2. Your FASTA file is indexed using `samtools`:

```bash
samtools faidx Data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna
```

3. Place the input files in the `Data/` directory:
   - BAM file: `SRR413984.sorted.NC_000001.10.bam`
   - BAM index: `SRR413984.sorted.NC_000001.10.bam.bai`
   - Reference FASTA: `GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna`

---

## Running the Script

Run the script with:

```bash
python bam_read.py
```

You should see a progress bar while start/end clusters are being matched. The output will be saved to:

```
detected_microdna_circles.txt
```

---

## Input Requirements

- **BAM file**: Contains aligned sequencing reads. Must be sorted and indexed.
- **FASTA file**: Reference genome for the target chromosome. Must be indexed.

---

## Output Details

The output file `detected_microdna_circles.txt` contains the following columns:

- **Chrom**: Chromosome name.
- **Start**: Start position of the detected circular DNA.
- **End**: End position of the detected circular DNA.
- **JunctionCount**: Total number of supporting soft-clipped reads.
- **Homology**: Microhomology overlap length between start and end sequences.
- **CompositeScore**: A score combining junction count and homology strength.

---

## Scoring Method

Each detected circle candidate is scored based on:
- The number of soft-clipped reads (`JunctionCount`).
- The level of **microhomology** between the start and end consensus sequences.

### Formula:

```python
CompositeScore = (Homology^2) × (JunctionCount^0.5)
```

This scoring scheme emphasizes **homology strength** while accounting for the number of supporting reads.

Only candidates with:
- **Microhomology ≥ 3 bp**
- **At least 10 soft-clipped reads total**

...are included in the final report.

---

## IGV Validation

To validate the results, open the output regions in **IGV** alongside the reference FASTA and BAM file. Look for:

- **Soft-clipped bases** (e.g., `2S40M` in the CIGAR string).
- **Opposite-direction clusters**.
- **Coverage resuming downstream** (suggests circular re-entry).

Recommended IGV settings:
- View → Color by → Strand
- View → Show soft-clipped bases
- View → Show all bases

---

## Output Example

```
Chrom           Start       End         JunctionCount   Homology    CompositeScore
NC_000001.10    24122430    24122550    890             3           268.50
```

---

## Script Functionality

The script performs the following steps:

1. **Soft-clip clustering**: Groups soft-clipped reads into start and end clusters based on their positions.
2. **Consensus sequence generation**: Derives consensus sequences for each cluster.
3. **Microhomology detection**: Aligns start and end consensus sequences to detect microhomology.
4. **Scoring**: Calculates a composite score for each candidate circle.
5. **Filtering**: Outputs only high-confidence candidates based on homology and junction count thresholds.

---

## Author Notes

This script was developed in a WSL Linux environment as part of a *Computational Genomics* course project to identify microDNA in real NGS data and analyze structural evidence for circular DNA formation.

---

## References

- [pysam docs](https://pysam.readthedocs.io/en/latest/)
- [parasail](https://github.com/jeffdaily/parasail)
- [microDNA discovery (Science, 2012)](https://www.science.org/doi/10.1126/science.1213307)
