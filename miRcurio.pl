#!/bin/perl
use lib "/home/alex/Desktop/miRcurio";
use strict;
use warnings;
use Getopt::Std;
use get_hairpins;
use get_degradome2;

use constant General_error => 1;

# produce error message if input data is incorrect or missing
sub usage {
	my ($message, $error) = @_;
	$error = General_error unless (defined $error);
	my @inputs = (
	'miRcurio: miRNA prediction pipeline (Dec 2019 Release)',
	'Usage: miRcurio [options] -s sRNA.fa -d degradome.fa -o output_dir',
	'Required :',
	'  -s sRNA_file.fa 	Path to the sRNA file',
	'  -d degradome_file.fa	Path to the degradome file',
	'  -o output_dir	Name for output files',
	'  -g genome.fa		Path to reference genome',
	'Optional :',
	'  -m mismatches	Mismatches for sRNA to its star sequence',
	'  -max 		Maximum number of times an sRNA can map to the genome',
	'  -flanking		Number of bases to select up- and down-stream of sRNA-genome match',
	'  -t temperature	Set the temperature for folding calculations - default is 20',
	'  -ffe			Set the maximum free folding energy - default -20',
	'  -seed_mismatch	maximum mismatches between the sRNA seed and degradome target',
	'  -seed_length		Set the length of the sRNA seed - default 9',
	'  -target_length	Set the length for the complement to the seed - default 12',
	'  -total_mismatch	Maximum number of mismatches between the miRNA and target',
	'  -max_targets		Maximum number of targets for each sRNA - defaul 100',

	);
	printf STDERR "$_\n" for @inputs;
	printf STDERR "\# $message\n";
	exit;
}

# define variables into a hash
my %opts;
getopts('s:d:o:m:g:max:flanking:t', \%opts);
usage ('Input sRNA file not specified!') unless (defined $opts{s});
usage ('Input degradome file not specified!') unless (defined $opts{d});
usage ('output directory not specified!') unless (defined $opts{o});
usage ('Genome file not specified!') unless (defined $opts{g});
$opts{m} = 6 unless (defined $opts{m} && $opts{m} =~ /^\d+$/);
$opts{max} = 99 unless (defined $opts{max} && $opts{max} =~ /^\d+$/);
$opts{flanking} = 300 unless (defined $opts{flanking} && $opts{flanking} =~ /^\d+$/);
$opts{t} = 20 unless (defined $opts{t} && $opts{t} =~ /^\d+$/);
$opts{ffe} = -20 unless (defined $opts{ffe} && $opts{ffe} =~ /^-\d+$/);
$opts{seed_mismatch} = 1 unless (defined $opts{seed_mismatch} && $opts{seed_mismatch} =~ /^\d+$/);
$opts{seed_length} = 9 unless (defined $opts{seed_length} && $opts{seed_length} =~ /^\d+$/);
$opts{target_length} = 12 unless (defined $opts{target_length} && $opts{target_length} =~ /^\d+$/);
$opts{total_mismatch} = 4 unless (defined $opts{total_mismatch} && $opts{total_mismatch} =~ /^\d+$/);
$opts{max_targets} = 100 unless (defined $opts{max_targets} && $opts{max_targets} =~ /^\d+$/);

# create output directory
system (mkdir $opts{o});

# use the short sequence aligner patman to map sRNA to the genome with no mismatches
Hairpin_analysis::patman ($opts{g}, $opts{s}, 0, "$opts{o}/temp_hairpin_patman");

# filter out sequences which map multiple times to the genome
Hairpin_analysis::filter_multimappers ("$opts{o}/temp_hairpin_patman", "$opts{o}/temp_hairpin_remove_duplicates", $opts{max});

# convert patman file to bedfile
Hairpin_analysis::patman_to_bed ("$opts{o}/temp_hairpin_remove_duplicates", $opts{flanking}, "$opts{o}/temp_hairpin_deduplicate.bed");

# filter overlapping sequences and select the most abundant sRNA
Hairpin_analysis::filter_overlap ("$opts{o}/temp_hairpin_deduplicate.bed", "$opts{o}/temp_hairpin_filter_overlap.bed");

# convert bedfile to fasta file
Hairpin_analysis::bed_to_fasta ($opts{g}, "$opts{o}/temp_hairpin_filter_overlap.bed", "$opts{o}/temp_hairpin_deduplicate.fasta");

# check flanking nucleotides for a possible complement to the sRNA
Hairpin_analysis::find_sRNA_complement ("$opts{o}/temp_hairpin_deduplicate.fasta", $opts{flanking}, $opts{m}, "$opts{o}/temp_hairpin_candidate_hairpin.bed");

# convert bed to fasta, this is required as previous sequences included up and down-stream flanking regions
Hairpin_analysis::bed_to_fasta ($opts{g}, "$opts{o}/temp_hairpin_candidate_hairpin.bed", "$opts{o}/temp_hairpin_candidate.fa");

# fold candidate precursors to confirm that they can form a viable hairpin
system ("RNAfold -i $opts{o}/temp_hairpin_candidate.fa -T $opts{t} --noPS >$opts{o}/temp_hairpin_folded");

# check the folded precursor to confirm that a viable hairpin can be formed.
Hairpin_analysis::hairpin_folding ("$opts{o}/temp_hairpin_folded", "$opts{o}/TEMPtemp_hairpin_candidate_precursors.bed", "$opts{o}/TEMP_hairpin_sRNA.fa", $opts{m}, $opts{ffe});

# Convert bedfile from folded hairpins back to a .fa file format
Hairpin_analysis::bed_to_fasta ($opts{g}, "$opts{o}/TEMPtemp_hairpin_candidate_precursors.bed", "$opts{o}/TEMPtemp_hairpin_candidate_precursors.fa");
##### same to here.

# map sRNA to putative hairpin precursors using patman.
Hairpin_analysis::patman ("$opts{o}/TEMPtemp_hairpin_candidate_precursors.fa", $opts{s}, 0, "$opts{o}/TEMP_hairpin_sRNA.patman");

# Check to confirm that an miRNA* can map to where the sRNA originated from.
Hairpin_analysis::check_complement ("$opts{o}/TEMP_hairpin_sRNA.patman", "$opts{o}/TEMP_candidate_hairpin");
######
######
print "up to degradome analysis\n";

# match the seed of sRNA to degradome sequences
Degradome_analysis::match_seed ("$opts{o}/TEMP_candidate_hairpin.sRNA.fa", $opts{d}, $opts{seed_length}, $opts{target_length}, $opts{seed_mismatch}, "$opts{o}/temp_degradome");


# Take degradome sequences which may form targets and map to the genome. From the alignments identify the total area an miRNA may map to.
Degradome_analysis::target_region ("$opts{o}/temp_degradome_degradome_targets.fa", $opts{g}, "$opts{o}/temp_degradome_target_region");

# Map the sRNA to the region from which a degradome sequence is derived to determine if this may form a putative target.
Degradome_analysis::cleavage_site_to_sRNA ("$opts{o}/temp_degradome_degradome_to_sRNA.fa", "$opts{o}/temp_degradome_target_region", $opts{total_mismatch}, "$opts{o}/temp_alignment_target");

# A final check to confirm sequences with a degradome target have a putative hairpin
# Running this script returns a final list of hairpin sequences, sRNA and their mapped locations.
Degradome_analysis::get_intersect ("$opts{o}/TEMP_candidate_hairpin", "$opts{o}/temp_alignment_target", $opts{max_targets}, "$opts{o}/FINAL_output"); 


# an example command line for running this script.
#perl miRcurio.pl -s sRNA.fa -g genome.fa -o output_directory -d degradome
