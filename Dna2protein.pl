#!/usr/bin/perl
system(`mkdir protein_sequence`);
use File::Basename;
@data= glob ("*.fas");
foreach $data (@data) {
$data1 = basename $data;

open (IN,"<$data");
open (OUT,">temp");
$first_line = 0;
	while(<IN>){
		 s/\r\n//;
		chomp;
		if($_=~/^\>/){
			if ($first_line eq 0) {
				print OUT "$_\n";
			}else{
			print OUT "\n",$_,"\n";}
			$first_line ++;
		}
		elsif($_){
			$_=~s/ //g;
			print OUT $_;
		}
		
	}

 close IN;
 close OUT;


open (IN,"<temp");
open (OUT,">./protein_sequence/$data1");
while (<IN>) {
	if ($_=~/\>/) {
		print OUT $_;
	}
	elsif ($_){
		 ($cds_seq) =$_=~/(.*)/;
         %standard_codon = ( 
 
GCT => 'A',
GCC => 'A',
GCA => 'A',
GCG => 'A',
GGT => 'G',
GGC => 'G',
GGA => 'G',
GGG => 'G',
GTT => 'V',
GTC => 'V',
GTA => 'V',
GTG => 'V',
GAT => 'D',
GAC => 'D',
GAA => 'E',
GAG => 'E',
CGT => 'R',
CGC => 'R',
CGA => 'R',
CGG => 'R',
CAA => 'Q',
CAG => 'Q',
CAT => 'H',
CAC => 'H',
CTT => 'L',
CTC => 'L',
CTA => 'L',
CTG => 'L',
CCT => 'P',
CCC => 'P',
CCA => 'P',
CCG => 'P',
AGA => 'R',
AGG => 'R',
AGT => 'S',
AGC => 'S',
AAT => 'N',
AAC => 'N',
AAA => 'K',
AAG => 'K',
ACT => 'T',
ACC => 'T',
ACA => 'T',
ACG => 'T',
ATT => 'I',
ATC => 'I',
ATA => 'I',
ATG => 'M',
TGT => 'C',
TGC => 'C',
TGA => '*',
TGG => 'W',
TTA => 'L',
TTG => 'L',
TTT => 'F',
TTC => 'F',
TCT => 'S',
TCC => 'S',
TCA => 'S',
TCG => 'S',
TAT => 'Y',
TAC => 'Y',
TAG => '*',
TAA => '*',
"---" => '-');
 $correspond_protein_seq="";

for (  $i=0; $i+3<=length($cds_seq); $i+=3 ) {        
     $codon = substr $cds_seq, $i, 3;    
	if ($standard_codon{ $codon } ne "") {
		 $aa= $standard_codon{ $codon };
		#print "$aa\t";
	}
	else {
         $aa= "X";
		#print "$aa\t";
	}
     print "$aa\t";
    $correspond_protein_seq .= $aa;
	#print  "$correspond_protein_seq"; 
}

print  OUT "$correspond_protein_seq\n"; 

	}

}
close IN;
close OUT;


}
