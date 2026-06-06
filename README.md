## (1) Extract_gff_features.py Usage Instructions
### Function: One-click output of 6 BED files (gene, exon, CDS, intron, gene upstream, distal intergenic) from the genome annotation file (gff) downloaded from NCBI.
## Usage command
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

## (2) to be continued
...
