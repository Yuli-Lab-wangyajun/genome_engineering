import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

def parse_gff(gff_file):
    gene_transcripts = defaultdict(dict)
    current_gene_id = None
    current_transcript_id = None

    with open(gff_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            seqid, source, feature, start, end, score, strand, phase, attributes = fields
            start, end = int(start), int(end)
            attr_dict = {}
            for attr in attributes.split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attr_dict[key.strip()] = value.strip()

            if feature == 'gene':
                if 'ID' in attr_dict:
                    current_gene_id = attr_dict['ID']
                else:
                    current_gene_id = None

            elif feature == 'mRNA':
                if 'ID' in attr_dict and current_gene_id:
                    current_transcript_id = attr_dict['ID']
                    gene_transcripts[current_gene_id][current_transcript_id] = {
                        'seqid': seqid,
                        'cds_length': 0,
                        'strand': strand,
                        'exons': []
                    }

            elif feature in ['exon', 'CDS']:
                if 'Parent' in attr_dict and current_gene_id:
                    parent_id = attr_dict['Parent']
                    if parent_id in gene_transcripts[current_gene_id]:
                        if feature == 'exon':
                            gene_transcripts[current_gene_id][parent_id]['exons'].append((start, end))
                        elif feature == 'CDS':
                            gene_transcripts[current_gene_id][parent_id]['cds_length'] += (end - start + 1)
    return gene_transcripts

def get_longest_transcript(gene_transcripts):
    longest_transcripts = {}
    for gene_id, transcripts in gene_transcripts.items():
        max_length = -1
        selected_tx = None
        for tx_id, tx_info in transcripts.items():
            if tx_info['cds_length'] > max_length:
                max_length = tx_info['cds_length']
                selected_tx = tx_id
        longest_transcripts[gene_id] = selected_tx
    return longest_transcripts

def extract_exon_sequences(genome, seqid, exons, strand):
    sequence = []
    sorted_exons = sorted(exons, key=lambda x: x[0], reverse=(strand == '-'))
    for start, end in sorted_exons:
        exon_seq = genome[seqid].seq[start-1:end]  # 提取 Seq 对象
        if strand == '-':
            exon_seq = exon_seq.reverse_complement()  # 直接操作 Seq 对象
        sequence.append(str(exon_seq))  # 转换为字符串
    return ''.join(sequence)

def main():
    parser = argparse.ArgumentParser(description='Extract longest transcript sequences from GFF and genome FASTA.')
    parser.add_argument('-g', '--gff', required=True, help='Input GFF file')
    parser.add_argument('-f', '--fasta', required=True, help='Genome FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    args = parser.parse_args()

    gene_transcripts = parse_gff(args.gff)
    print(f"[DEBUG] Parsed {len(gene_transcripts)} genes")

    longest_transcripts = get_longest_transcript(gene_transcripts)
    print(f"[DEBUG] Found {len(longest_transcripts)} longest transcripts")

    genome = SeqIO.index(args.fasta, 'fasta')

    with open(args.output, 'w') as f_out:
        for gene_id, tx_id in longest_transcripts.items():
            if not tx_id:
                continue
            tx_info = gene_transcripts[gene_id][tx_id]
            seqid = tx_info['seqid']
            if seqid not in genome:
                print(f"[WARNING] Missing seqid: {seqid}")
                continue

            full_sequence = extract_exon_sequences(genome, seqid, tx_info['exons'], tx_info['strand'])
            header = f">{gene_id}"
            f_out.write(f"{header}\n{full_sequence}\n")
    print(f"[INFO] Output saved to {args.output}")

if __name__ == '__main__':
    main()