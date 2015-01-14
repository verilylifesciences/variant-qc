#!/usr/bin/perl
#
# This script takes a vcf file as input and outputs a new vcf file with only positions with transition or transversion mutations.  The header remains intact.
#
# Usage: pull_titv.pl [vcf_file] 

use strict;
use warnings;

# VCF file from command line
my $file = $ARGV[0];

# Define transition and transversion mutations
my %transitions = ( 
	A => 'G',
	G => 'A',
	C => 'T',
	T => 'C',
);

my %transversions = (
	A => ['C','T'],
	C => ['A','G'],
	T => ['G','A'],
	G => ['T','C'],
);

# Open input file and check variants
open(DATA,$file);

for(<DATA>) {
	my $line = $_;

	if ( $line =~ /^#/ ) {
		print $line;
		next;
	}

	my @row = split(/\s+/,$line);
	my $ref = $row[3];
	my $alt = $row[4];
	
	my $match = 0;

	if ( exists $transitions{$ref}  && $transitions{$ref} eq $alt ) {
		$match = 1;
	}
	
	elsif ( exists $transversions{$ref} ) {
		
		my %alts = map { $_ => 1 } @{$transversions{$ref}};
		if ( exists $alts{$alt} ) {
			$match = 1;
		}
	}
	
	if ( $match == 1 ) {
		print $line;
	}
}
