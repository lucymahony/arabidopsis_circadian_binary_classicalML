import pandas as pd
from collections import Counter
import sys
from pyfaidx import Fasta
from itertools import product
# class 1 (circadian) from class 0 (non-circadian).
# “TATTGC” for counts from the promoter and “TATTGC.1” for counts from the mRNA. 


def count_kmers(sequence, k):
    """Count the occurrences of k-mers in a given sequence."""
    sequence = sequence.upper() # Convert to uppercase to make the count case-insensitive
    valid_bases = {'A', 'C', 'T', 'G'}
    return Counter(
        kmer for i in range(len(sequence) - k + 1)
        if set(kmer := sequence[i:i+k]).issubset(valid_bases))


def generate_matrix(genes_labels_file, k):
    # Generate a DataFrame with transcript name, the circadian label, and the k-mer count columns currently set to 0
    genes_labels_df = pd.read_csv(genes_labels_file, sep=',')
    genes_labels_df = genes_labels_df[['Transcript', 'circadian']]
    genes_labels_df = genes_labels_df.rename(columns={'Transcript': 'transcript', 'circadian': 'label'})
    bases = ['A', 'C', 'G', 'T']
    kmers = [''.join(p) for p in product(bases, repeat=k)]
    all_kmer_cols = kmers + [f"{kmer}.1" for kmer in kmers]
    kmer_df = pd.DataFrame(0, index=genes_labels_df.index, columns=all_kmer_cols)
    genes_labels_df = pd.concat([genes_labels_df, kmer_df], axis=1)
    return genes_labels_df


def fill_in_kmers(genes_labels_df, transcripts_seqs, promoter_seqs, k):
    # Fill in the k-mer counts for each transcript and promoter sequence
    # Convert to uppercase to make the search case-insensitive
    # if the transcript isnt found, try appending ".1" to the transcript ID to search for that transcript. 
    for index, row in genes_labels_df.iterrows():
        original_id = row['transcript']
        transcript_id = original_id.upper()  # Convert to uppercase to make the search case-insensitive

        # Attempt with original ID
        found_transcript = transcript_id in transcripts_seqs
        found_promoter = transcript_id in promoter_seqs

        # If not found, try appending ".1"
        if not found_transcript and not found_promoter:
            transcript_id_with_suffix = f"{transcript_id}.1"
            found_transcript = transcript_id_with_suffix in transcripts_seqs
            found_promoter = transcript_id_with_suffix in promoter_seqs
            if found_transcript or found_promoter:
                transcript_id = transcript_id_with_suffix

        if found_transcript:
            transcript_seq = str(transcripts_seqs[transcript_id])
            transcript_kmer_counts = count_kmers(transcript_seq, k)
            for kmer in transcript_kmer_counts:
                genes_labels_df.at[index, kmer] = transcript_kmer_counts[kmer]
        else:
            print(f"Transcript ID {original_id} not found in transcript sequences (even with .1).")

        if found_promoter:
            promoter_seq = str(promoter_seqs[transcript_id])
            promoter_kmer_counts = count_kmers(promoter_seq, k)
            for kmer in promoter_kmer_counts:
                genes_labels_df.at[index, f"{kmer}.1"] = promoter_kmer_counts[kmer]
        else:
            print(f"Promoter ID {original_id} not found in promoter sequences (even with .1).")

    return genes_labels_df


if __name__ == "__main__":
    # Input file paths - test files
    # genes_labels_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/test_pnas.csv"     
    # promoter_fasta =  "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/test_promoters.fa" 
    # transcript_fasta =  "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/test_transcripts.fa"
    
    # Input file paths
    genes_labels_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/raw/pnas.2103070118.sd01.csv"     
    promoter_fasta =  "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/promoters.fa" 
    transcript_fasta =  "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/transcripts.fa"

    # Output file paths
    matrix_output_file = "/ei/projects/c/c3109f4b-0db1-43ec-8cb5-df48d8ea89d0/scratch/repos/arabidopsis_circadian_binary_classicalML/data/processed/promoters_transcripts_6_mers_labelled.csv" # Column names transcript, kmer ,label  

    # Read in parameters
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 6 # default is 6 

    # Read in promoter fasta, needs to be cleaned as :: then info is added to the fasta header 
    promoter_seqs = Fasta(promoter_fasta)
    promoter_seqs_cleaned = {}
    for full_id in promoter_seqs.keys():
        clean_id = full_id.split("::")[0]
        promoter_seqs_cleaned[clean_id] = str(promoter_seqs[full_id])
    promoter_seqs = promoter_seqs_cleaned

    transcript_seqs = Fasta(transcript_fasta)

    genes_labels_df = generate_matrix(genes_labels_file, k)
    print('Filling in k-mers ... ')
    genes_labels_df = fill_in_kmers(genes_labels_df, transcript_seqs, promoter_seqs, k)
    genes_labels_df.to_csv(matrix_output_file, index=False)
    print(f"Matrix saved to {matrix_output_file}")