#!/bin/perl
use lib "/home/alex/Desktop/miRcurio";
use strict;
use warnings;
use Getopt::Std;
use get_hairpins;
use get_degradome;

my ($input) = @ARGV;

# script to process cleavage_site_to_sRNA in parallel. 
# You will need to modify input sequences so that filepaths are specific to your computer.
# This is run using gnu parallel 
Degradome_analysis::cleavage_site_to_sRNA ("known_miRNA2/temp_degradome_degradome_to_sRNA.fa", $input, 4, "known_miRNA2/temp_alignment_target");

# This is the command line input that is used to run this script in parallel
#cat output_wheat/temp_degradome_target_region | parallel -k -j 100 --pipe -N50000 perl run_parallel.pl -
