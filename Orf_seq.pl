#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;

#set parameter
my $version=1.00;
my %opts;
GetOptions(\%opts,"r=s","w=s","h");
if (!(defined $opts{w} && defined $opts{r}) || defined $opts{h})
{
        &usage;
}
open Input,"$opts{r}";   #fas_l.list
open Output,">$opts{w}.c.fas";
open Check,">check";
my $sp="";
my $seq="";
my $len=0;
my $check=0;
while (<Input>){
	s/\n|\r//g;
  if(/>(\w+)_\d+ \[(\d+) - (\d+)\]/) {
  	print Check "$1\t$2\t$3\t";
  	if (($1 ne $sp) and ($sp ne "")){
  		print Output ">$sp\n$seq\n";
  		$seq="";
  		$check=1;
        if ($3>$2) {
            $len=$3-$2;
        }
        else {
            $len=$2-$3;
        }
  		
  	}
  	else {
		if ($3>$2) {
  			if (($3-$2)>$len) {
  				$check=1;
  				$seq="";
				$len=$3-$2;
  			}
  			else {$check=0;}  
		}
		else {
 			if (($2-$3)>$len) {
  				$check=1;
  				$seq="";
				$len=$2-$3;
  			}
  			else {$check=0;}  	
		}		
  	}
  	$sp=$1;
  	print Check "$len";
  }
  else {
  	if ($check==1){$seq.=$_;}
  }
}
print Output ">$sp\n$seq\n";
close Input;
close Output;

