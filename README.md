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

## (3) Obtain the ORF from the longest transcript file (Dna2protein.pl and Orf_seq.pl)
```
# Change to the working directory
cd path/test

# 1. Rename sequence headers: replace hyphens with underscores and create a name mapping file
sed -i.bak '/^>/ s/-/_/g' human_longestTRAN.fa;grep '^>' human_longestTRAN.fa.bak | sed 's/^>//' > old_names.txt;grep '^>' human_longestTRAN.fa | sed 's/^>//' > new_names.txt;paste old_names.txt new_names.txt > Name_align.txt;rm old_names.txt new_names.txt human_longestTRAN.fa.bak

# 2. Create a directory to store individual sequence files. Split the multi‑FASTA file into separate files by sequence name (using UCSC tools)
mkdir split
faSplit byname human_longestTRAN.fa split/

# 3. For each input sequence, predict all ORFs (≥120 aa, both strands) using EMBOSS getorf
raw_file_list=`ls split/`
mkdir orf
for file in $raw_file_list
do
new_file=${file%.*}
/work/home/acm2g39zj0/software/EMBOSS-6.5.7/emboss/getorf -sequence split/$file orf/$new_file.orf.fas -find 2 -minsize 120 -reverse
done

# 4. For each ORF file, extract the single longest ORF (using an external function and a Perl script)
mkdir orf_longest
orf_file_list=`ls orf/`
for orf_file in $orf_file_list
do
get_the_longest_orf
new_orf_file=${orf_file%.*}
perl $pl_path/Orf_seq.pl -r orf/$orf_file -w orf_longest/$new_orf_file
done

# 5. Translate the nucleotide ORFs into protein sequences; 
cd orf_longest
perl $pl_path/Dna2protein.pl

# 6. Go back, concatenate all protein FASTA files, remove any lines that are just '>', and save as final file
cd ..
cat orf_longest/protein_sequence/*.fas | sed '/^>$/d' > human_longestTRAN_protein.orf.all.fa

```

## (4) to be continued
...
