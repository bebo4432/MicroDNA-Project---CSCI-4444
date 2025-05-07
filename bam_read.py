import pysam
from collections import Counter
import numpy as np
import parasail
from collections import defaultdict
from tqdm import tqdm


def find_microhomology(seq1, seq2, min_h=3, max_h=10):
    alignment = parasail.sw_trace(seq1, seq2, 10, 1, parasail.dnafull)
    score = alignment.score
    ref_aligned = alignment.traceback.ref
    query_aligned = alignment.traceback.query
    matches = sum(1 for r, q in zip(ref_aligned, query_aligned) if r == q and r != '-')
    if min_h <= matches <= max_h:
        return matches
    return 0

def fetch_reference_sequence(ref_path, chrom, start, end):
    if start >= end:
        return ""
    try:
        ref = pysam.FastaFile(ref_path)
        return ref.fetch(chrom, start, end)
    except Exception as e:
        print(f"[!] Failed to fetch {chrom}:{start}-{end}: {e}")
        return ""


def count_matching_chars(char, char_list):
    return sum(1 for c in char_list if c == char)

def most_common_char(char_list):
    return Counter(char_list).most_common(1)[0][0]

def get_consensus(seq_list, min_chars=15, min_density=0.75):
    consensus = {}
    S = ""
    for seq in seq_list:
        for i in range(len(seq)):
            if i not in consensus:
                consensus[i] = []
            consensus[i].append(seq[i])
    for i in sorted(consensus.keys()):
        mcc = most_common_char(consensus[i])
        mcc_count = count_matching_chars(mcc, consensus[i])
        density = mcc_count / len(consensus[i])
        if mcc_count >= min_chars and density >= min_density:
            S += mcc
    return S


def get_clipped_seq(read):
    clipped_seq = ""
    pos = 0
    for op, length in read.cigartuples:
        if op == 4:  # soft clip
            clipped_seq += read.query_sequence[pos:pos+length]
            pos += length
        elif op == 0:  # match
            pos += length
    return clipped_seq

# Cluster settings
ref_fasta = "Data/GCF_000001405.13_GRCh37_genomic.NC_000001.10.fna"
bam_file = "Data/SRR413984.sorted.NC_000001.10.bam"
min_cluster_support = 5  # minimum number of clipped reads to consider
output_file = "detected_microdna_circles.txt"

# Cluster soft-clips
start_clusters = defaultdict(list)
end_clusters = defaultdict(list)

with pysam.AlignmentFile(bam_file, "rb") as bam:
    print("BAM references:", bam.references)
    for read in bam:
        if not read.cigartuples:
            continue

        chrom = read.reference_name
        if read.cigartuples[0][0] == 4:  # softclip at start
            pos = read.reference_start
            start_clusters[(chrom, pos)].append(get_clipped_seq(read))
        elif read.cigartuples[-1][0] == 4:  # softclip at end
            pos = read.reference_end
            end_clusters[(chrom, pos)].append(get_clipped_seq(read))


# estimate total comparisons
total = len(start_clusters) * len(end_clusters)

with open(output_file, "w") as out:

    out.write("Chrom\tStart\tEnd\tJunctionCount\tHomology\tCompositeScore\n")

    for (chrom, s_pos), s_clips in tqdm(start_clusters.items(), desc="Start clusters"):
        for (chrom2, e_pos), e_clips in end_clusters.items():
            if chrom != chrom2:
                continue
            if len(s_clips) < min_cluster_support or len(e_clips) < min_cluster_support:
                continue
            if abs(e_pos - s_pos) > 2000:
                continue

            start_consensus = get_consensus(s_clips)
            end_consensus = get_consensus(e_clips)

            if not start_consensus or not end_consensus:
                continue

            ref_seq = fetch_reference_sequence(ref_fasta, chrom, s_pos, e_pos)
            if not ref_seq:
                continue
            homology = find_microhomology(start_consensus[::-1], end_consensus)

            num_junctions = len(s_clips) + len(e_clips)

            # Optional: Only report good circles
            if homology >= 3 and num_junctions >= 10:
                composite_score = num_junctions * homology
                out.write(f"{chrom}\t{s_pos}\t{e_pos}\t{num_junctions}\t{homology}\t{composite_score}\n")

                
                

if __name__ == "__main__":
    print("Starting microDNA detection...")
    # Code above runs automatically
    print("Finished! Results in:", output_file)








