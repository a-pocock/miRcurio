#!/bin/perl
#
package Hairpin_analysis;
use strict;
use warnings;
use File::Temp qw(tempfile);

# run patman, a short sequence aligner
# gaps are allowed, however, these cannot exceed mismatches
sub patman {
	my ($database, $pattern, $mismatch, $patout) = @_;
	system ("patman -D $database -P $pattern -e $mismatch -o $patout -g 2");
	return;
}

# convert patman alignment to bed format
# take flanking regions up and downstream of alignment
# information about the alignment and orientation is also retained
sub patman_to_bed {
	my ($patman_file, $flanking, $bed_out) = @_;
	open (IN, $patman_file) or die "can't open $!";
	open (OUT, ">>$bed_out") or die "can't open $!";
	# open patman file, read each line and convert to .bed.
	foreach (<IN>) {
		$_ =~ /(\S+)\t(\S+)\t(\d+)\t(\d+)\t(.)\t(\d+)/;
		my $left_position = $3-1;
		my $left_start = $left_position - $flanking; my $right_stop = $4 + $flanking;
		# check values are not negative and can form true bed values
		if ($left_start < 0) {
			$left_start = 0;
		}
		my $bedfile = "$1\t$left_start\t$right_stop\t$1_$2_$left_position\_$4_$5_$6\n";
		print OUT $bedfile;
	}
}

# read through bedfile and convert to fasta using getfasta from bedtools
sub bed_to_fasta {
	my ($genome, $bedfile, $out_fasta) = @_;
	system "bedtools getfasta -fi $genome -bed $bedfile -s -name -fo $out_fasta";
}

# take region around mapped sRNA and check to see if a  reverse complement to the sRNA maps either up or downstream from the miRNA. 
# if there is a reverse complement to the miRNA, output a bedfile with coordinates including the complement to the miRNA
sub find_sRNA_complement {
	my ($fasta_file, $offset, $mismatches, $output_bed) = @_;
	open (IN, $fasta_file) or die "can't open $!";
	open (BEDOUT, ">>$output_bed") or die "can't open $!";
	while (<IN>) {
		if ($_ =~ /^>\(([^)]+)\)_(\w+)_(\d+)_(\d+)_(.)_(.)/) {
			# get values for patman comparison
			my $chr=$1;my $sRNA=$2;my $start_position=$3;my $stop_position=$4;my $orient=$5;
			my $header = "($1)_$2_$3_$4_$5_$6";
			my $seq = (<IN>);
			my $sRNA_length = length($sRNA);
			# create temporary directories for patman input and output databases
			my $tempout = File::Temp->new();
		        my $temp_database = File::Temp->new();
        		my $temp_sRNA = File::Temp->new();
			my $lowstart = $start_position; my $highstop = $stop_position;
			open (PATTERN, ">>$temp_sRNA") or die "can't open $!";
			open (DATABASE, ">>$temp_database") or die "can't open $!";
			# populate temporary files for make database and pattern for patman
			print DATABASE ">$sRNA\n$seq"; print PATTERN ">$sRNA\n$sRNA";
			patman ($temp_database, $temp_sRNA, $mismatches, $tempout);
			open (PAT_OUT, $tempout) or die "can't open $!";
			while (<PAT_OUT>) {
				if ($_ =~ /(\S+)\t(\S+)\t(\S+)\t(\S+)\t(.)/) {
					# check alignments and ensure sRNA is in opposite orientation to its complement
					my $alignment_start = $3-1;
					if ($5 ne $orient) {
						if ($4 < $offset) {
							my $start = $start_position-($offset-$alignment_start);
							if ($start < $lowstart) {
								$lowstart = $start;
							}
						}
						if ($alignment_start > ($offset+$sRNA_length)) {
							my $end = ($4-$offset) + $start_position;
							if ($end > $highstop) {
								$highstop = $end;
							}
						}
					}
				}
			}
			if ($lowstart != $start_position) {
				print BEDOUT "$chr\t$lowstart\t$stop_position\t$header\t1\t$orient\n";
			}
			if ($highstop != $stop_position) {
				print BEDOUT "$chr\t$start_position\t$highstop\t$header\t1\t$orient\n";
			}
		}
	}
}

# remove sequences which map greater than X number of times to the genome (default 100)
sub filter_multimappers {
	my ($patman_file, $output, $maximum_count) = @_;
	open (IN, $patman_file) or die "can't open $!";
	open (OUT, ">>$output") or die "can't open $!";
	# create hash and generate a count for each sRNA
	my %sequences;
	while (<IN>) {
		if ($_ =~ /^\S+\t([^-]+)/) {
			$sequences{$1}++;
		}
	}
	# iterate over the patman input file again and print values at or below the maximum count
	open (IN, $patman_file) or die "can't open $!";
	while (<IN>) {
		if ($_ =~ /^\S+\t([^-]+)/) {
			if ($sequences{$1} <= $maximum_count) {
				print OUT "$_";
			}
		}
	}
}

# filter mapped sequences for overlapping regions in the genome. 
# Where multiple sequences map to the genome the most abundant sRNA is retained for further analysis
sub filter_overlap {
	my ($input_patman, $output_patman) = @_;
	# begin by sorting sequences on genomic location
	system "sort -k 1,1 -k 2,2n -k 3,3n $input_patman > tmp_sorted_matches";
	open (IN, "tmp_sorted_matches") or die "can't open $!";
	my @alignments;
	# set initial values for sRNA count and the end position.
	my $count = 0; my $end = 0; my $chr = "chr"; my $start;
	open (OUT, ">>$output_patman") or die "can't open $!";
	# read through each line of the input file, check if there is overlap between sequences
	# if overlap determine most abundant sequence using count
	while (<IN>) {
		$_ =~ /^(\S+)\t(\S+)\t(\S+)\t\1(_[A-Z]+)-(\d+)(.+)/;
		# check if using the same chromosome and if there is overlap
		if (($1 eq $chr) and ($2 < $end)) {
			$end = $3;
			if ($5 > $count) {
				# replace current data in @alignments
				$chr = $1; $start = $2; $end = $3; $count = $5;
				@alignments = ("$1\t$2\t$3\t($1)$4$6");
			}
			elsif ($5 == $count) {
				# append to @alignments
				$chr = $1; $start = $2; $end = $3; $count = $5;
				unshift @alignments, "$1\t$2\t$3\t($1)$4$6";
			}
		}
		else {
			# print each item in @alignments
			foreach (@alignments) {
				print OUT "$_\n";
			}
			# create new list called @lines
			$chr = $1; $start = $2; $end = $3; $count = $5;
			@alignments = ("$1\t$2\t$3\t($1)$4$6");
		}
	}
	# at end of list ensure final value of @alignments is printed
	foreach (@alignments) {
		print OUT "$_\n";
	}
}

# take hairpins folded using RNAfold and check that a viable hairpin can be made
sub hairpin_folding {
	my ($RNAfold, $output, $output_fasta, $mismatch, $minimum_ffe) = @_;
	open (IN, $RNAfold) or die "can't open $!";
	open (OUT, ">>$output") or die "can't open $!";
	open (OUT_sRNA, ">>$output_fasta") or die "can't open $!";
	# read each header, extract key data and then read hairpin sequence and folded structure
	while (<IN>) {
		if ($_ =~ /^>\(([^)]+)\)_(\w+)_(\d+)_(\d+)_(.)_\d+::\1:(\d+)-(\d+)/){
			my $chromosome = $1; my $sRNA = $2; my $sRNA_start = $3; my $sRNA_end = $4; my $orientation = $5; my $hairpin_start = $6; my $hairpin_end = $7; my $hairpin_orientation;
			my $hairpin_sequence = <IN>;
			my $fold = <IN>; chomp $fold;
			# get the structure and folding energy
			# if folding energy is above the threshold skip current sRNA and go onto the next
			unless ($fold =~ /(\S+) \( *(-\d+\.\d+)\)/) {
				next;
			}
			my @hairpin = ($1, $2);
			my $result;
			if (int($hairpin[1]) >= $minimum_ffe) {
				next;
			}
			# determine if miRNA is on the 5' or 3' end of the precursor. 
			# if reading the precursor from the left the sRNA maps to the end
			# reverse the hairpin sequence before running checkhairpin
			if ($sRNA_start == $hairpin_start) {
				if ($orientation eq '-') {
					my $revhairpin = reverse_hairpin ($hairpin[0]);
					$result = check_hairpin ($sRNA, $revhairpin, $mismatch);
				}
				elsif ($orientation eq '+') {
					$result = check_hairpin ($sRNA, $hairpin[0], $mismatch);
				}
			}
			elsif ($sRNA_end == $hairpin_end) {
				if ($orientation eq '-') {
					$result = check_hairpin ($sRNA, $hairpin[0], $mismatch);
				}
				elsif ($orientation eq '+') {
					my $revhairpin = reverse_hairpin ($hairpin[0]);
					$result = check_hairpin ($sRNA, $revhairpin, $mismatch);
				}
			}
			# if the output of checkhairpin is true print key values from the input line
			if ($result eq "true") {
				print OUT "$chromosome\t$hairpin_start\t$hairpin_end\t$sRNA\_($chromosome)\_$sRNA_start\_$sRNA_end\_$hairpin_start\_$hairpin_end\_$orientation\t1\t$orientation\n";
				print OUT_sRNA ">$sRNA\n$sRNA\n";
				#				print OUT "$sRNA\t$chromosome\t$sRNA_start\t$sRNA_end\t$hairpin_start\t$hairpin_end\t$orientation\n"
			}
		}
	}
}

# takes folded structures from RNAfold and reverses the sequence.
# This includes reversing the direction of brackets contained within the secondary structure
sub reverse_hairpin {
	my ($folded_sequence) = @_;
	my $new_sequence = reverse ($folded_sequence);
	$new_sequence =~ tr/()/)(/;
	return $new_sequence;
}

# take folded hairpin and check that the precursor forms a viable hairpin structure
sub check_hairpin {
	my ($miRNA, $folded_hairpin, $mismatches) = @_;
	my $count = 0; my $miR_match = 0; my $star_count =0;
	my @hairpinstructure = split('', $folded_hairpin, length($folded_hairpin));
	my $miRNA_length = length($miRNA);
	my $hairpin_length = length($folded_hairpin);
	# check the miRNA does not fold back on itself
	while ($count < $miRNA_length) {
		if ($hairpinstructure[$count] eq "(") {
			$miR_match++;
		}
		elsif ($hairpinstructure[$count] eq ")") {
			return 'false';
		}
		$count++;
	}
	# check that the miRNA-star does not fold back on itself
	$count = 0;
	while ($star_count < $miR_match) {
		$count++;
		my $star_position = $hairpin_length-$count;
		if ($hairpinstructure[$star_position] eq ")") {
			$star_count++;
		}
		elsif ($hairpinstructure[$star_position] eq "(") {
			return 'false';
		}
	}
	# check the miRNA and star match to each other and have less than $mismatches between them
	if (($count >= $miRNA_length-$mismatches) && ($count <= $miRNA_length+$mismatches) && ($miR_match >= $miRNA_length-$mismatches)) {
		my $position = $miRNA_length;
		my $hairpin_value = 0;
		while ($position < $hairpin_length-$count) {
		#		while ($position < $hairpin_length-$star_count) {
			if ($hairpin_value < 0) {
				return "false";
			}
			elsif ($hairpinstructure[$position] eq '(') {
				$hairpin_value++;
			}
			elsif ($hairpinstructure[$position] eq ')') {
				$hairpin_value--;
			}
			$position++;
		}
		if ($hairpin_value == 0) {
			return "true";
		}
		else {
			return "false";
		}
	}
	else {
		return "false";
	}
}

# take patman alignment for all sRNA against predicted precursors if there is an sRNA with a perfect alignment then output the sequence as a potential precursor.
sub check_complement {
	my ($input_patman, $output) = @_;
	open (IN, "$input_patman") or die "can't open $!";
	my $chromosome; my $lowstart = 20; my $highstop = 15; my $seq = "ATCG"; my $start_site; my %hash; my %miRNA_seq; my $lowend; my $highstart;
	# check if alignment is in the same precursor or a new sequence
	# check that there is an sRNA match at either end of the putative hairpin
	while (<IN>) {
		$_ =~ /^([^_]+)([^:]+)::([^:]+):(\d+)-(\d+)\S+\t\w+-\d+\t(\d+)\t(\d+)/;
		my $length = $5-$4;
		my $header = "$1$2";
		if ($seq eq $1 && $chromosome eq $3 && $start_site == $4) {
			if ( $lowstart > $6) {$lowstart = $6; $lowend = $7;}
			if ( $highstop < $7) {$highstop = $7; $highstart = $6;}
			if (($lowstart < 15) and ($highstop > $length-15) and ($lowend < $highstart)) {
				if ($seq =~ /A/ && $seq =~ /T/ && $seq =~ /C/ && $seq =~ /G/) {
					$hash{$header} = $seq;
					$miRNA_seq{$seq} = $seq;
				}
			}
		}
		else {
			$seq = $1; $chromosome = $3; $start_site = $4; $lowstart = $6; $highstop = $7; $lowend =  $7; $highstart = $6;
		}
	}
	open (OUT, ">>$output") or die "can't open $!";
	foreach my $key (keys %hash) {
		print OUT ">$key\n$hash{$key}\n";
	}
	open (OUT2, ">>$output.sRNA.fa") or die "can't open $!";
	foreach my $sRNA (keys %miRNA_seq) {
		print OUT2 ">$sRNA\n$sRNA\n"; 
	}
}

1;
