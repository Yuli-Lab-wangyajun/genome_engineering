## (1) Extract_gff_features.py Usage Instructions
### Function: One-click output of 6 BED files (gene, exon, CDS, intron, gene upstream, distal intergenic) from the genome annotation file (gff) downloaded from NCBI.
```
usage: Extract_gff_features.py [-h] [-o PREFIX] [--chromosome CHROMOSOME]
                               [-up UPSTREAM]
                               gff

Example: python Extract_gff_features.py input.gff -o out --up 2000
-o: Prefix for output file names
-up: Upstream flanking length, default 2000bp
--chromosome: Optional, extract only a single chromosome
```
## Output results
```
genome_gene.bed、genome_exon.bed、genome_CDS.bed、genome_intron.bed、genome_upstream2000.bed、genome_distal_intergenic.bed
```

## (2) Extract_longest_transcript.py Usage Instructions
### Function: This script extracts the longest transcript (by total CDS length) for each gene from a GFF annotation file and a genome FASTA file. It then outputs the concatenated exon sequences (correctly oriented according to strand) as a FASTA file.
```
usage: python extract_longest_transcript.py -g annotation.gff -f genome.fasta -o longest_transcripts.fasta

-g, --gff	Yes	Input GFF file (version 3) with gene, mRNA, exon, and CDS features
-f, --fasta	Yes	Genome FASTA file (indexed automatically)
-o, --output	Yes	Output FASTA file containing the longest transcript sequences
```

## (3) to be continued
...
