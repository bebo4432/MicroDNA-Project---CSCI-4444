# microDNA Detector

This project implements a Python-based algorithm to detect **microDNA** (extrachromosomal circular DNA) in a cancer cell line from aligned sequencing data. The algorithm identifies and scores candidate circular DNA junctions using soft-clipped reads from a BAM file and validates them using microhomology detection.

---

## 📂 Project Structure

```
├── bam_read.py                  # Main microDNA detection script
├── detected_microdna_circles.txt # Output: list of predicted circles
└── Data/
    ├── SRR413984.sorted.NC_000001.10.bam   # Input BAM file
    └── GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna  # Chromosome 1 FASTA
```

---

## 📦 Requirements

This script is intended to run in a **Python 3 WSL (Linux)** environment.

Install the required packages:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install pysam parasail tqdm numpy
sudo apt install samtools  # for indexing FASTA
```

---

## 🔧 Setup

Ensure the following:

1. Your BAM file is sorted and indexed (`.bam` + `.bai`)
2. Your FASTA file is indexed using `samtools`:

```bash
samtools faidx Data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna
```

---

## ▶️ Running the Script

```bash
python bam_read.py
```

You should see a progress bar while start/end clusters are being matched. Output will be saved to:

```
detected_microdna_circles.txt
```

Each line reports:
- Chromosome
- Start position
- End position
- Number of supporting soft-clipped reads (`JunctionCount`)
- Microhomology overlap length (`Homology`)
- Composite score = `JunctionCount × Homology`

---

## 📊 Scoring Method

For each candidate circle, the script calculates:

- **Consensus** start and end sequences using high-density agreement across reads
- **Microhomology** between start and end tags using Smith-Waterman alignment (`parasail`)
- **Score** = total number of start/end tags × homology length

Only circles with:
- Homology ≥ 3 bp
- At least 10 supporting soft-clipped reads

...are included in the final output.

---

## 🔬 Validation

Open the output positions in **IGV** using the input BAM and FASTA to visualize circular junctions. Look for:

- Soft-clipped reads at both ends
- Reversely aligned sequences matching circularization

---

## 📁 Output Example

```
Chrom           Start       End         JunctionCount   Homology    CompositeScore
NC_000001.10    121484880   121485434   6567            3           19701
```

---

## ✍️ Author Notes

This script was developed in a WSL Linux environment as part of the *Computational Genomics* course project to identify microDNA in sequencing data and evaluate circular DNA formation using real alignment evidence.

---

## 🧠 References

- [pysam docs](https://pysam.readthedocs.io/en/latest/)
- [parasail](https://github.com/jeffdaily/parasail)
- [microDNA discovery (Science, 2012)](https://www.science.org/doi/10.1126/science.1213307)