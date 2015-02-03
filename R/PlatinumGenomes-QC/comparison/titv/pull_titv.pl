#!/usr/bin/env perl
#
# Copyright 2014 Google Inc. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
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
