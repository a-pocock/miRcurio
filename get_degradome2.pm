#!/bin/perl

package Degradome_analysis;
use strict;
use warnings;
use File::Temp qw(tempfile);

# get the first X sequences of the sRNA and the first Y sequences of the degradome. Use the short read aligner patman to match these two datasets together.
# output degradome sequences which are targeted by sRNA
sub match_seed {
	my ($sRNA, $degradome, $seed_length, $deg_target_length, $mismatches, $output) = @_;
	open (SRNA_IN, $sRNA) or die "can't open $!";
	my $temp_sRNA_seed = File::Temp->new();
	open (SEED, ">>$temp_sRNA_seed") or die "can't open $!";
	while (<SRNA_IN>) {
		if ($_ !~ /^>/) { 
			if ($_ =~ /^(\w{$seed_length})/) {
				print SEED ">$_";
				print ">$_";
				my $seq = $1;
				$seq =~ tr/ATGC/TACG/; 
				$seq = reverse($seq);
				print SEED "$seq\n";
				print "$seq\n";
			}
		}
	}
	open (DEG, $degradome) or die "can't open $!";
	my $temp_deg_target = File::Temp->new();
	open (DEG_TARGET, ">>$temp_deg_target") or die "can't open $!";
	while (<DEG>) {
		if ($_ =~ /^>/) {
			print DEG_TARGET "$_";
		}
		elsif ($_ =~ /(\w{$deg_target_length})/) {
			print DEG_TARGET "$1\n";
		}
	}
	my $temp_alignment = File::Temp->new();

	### maybe use system so the -s option can be added.
	###
	####
	print "got here\n";
	system ("patman -D $temp_deg_target -P $temp_sRNA_seed -e $mismatches -o $temp_alignment -s");
	#	Hairpin_analysis::patman ($temp_deg_target, $temp_sRNA_seed, $mismatches, $temp_alignment);
	open (PATMAN, $temp_alignment) or die "can't open $!";
	my %targets;
	my %targets_to_sRNA;
	open (OUT, ">>$output\_degradome_targets.fa") or die "can't open $!";
	open (OUT2, ">>$output\_degradome_to_sRNA.fa") or die "can't open $!";
	my $counter;
	while (<PATMAN>) {
		$_ =~ /^([^-]+)-(\d+)\t(\w+)/;
		$targets{"$1-$2"} = $1;
		$targets_to_sRNA{"$1-$2\t$3"} = "$3";
		$counter ++;
		print "$counter\n";
	}
	foreach my $degradome_seq (keys %targets) {
		print OUT ">$degradome_seq\n$targets{$degradome_seq}\n";
	}	

	foreach my $sRNA_deg (keys %targets_to_sRNA) {
		print OUT2 "$sRNA_deg\n";
	}
}

# get regions surrounding the cleavage site of degraded sequences
#
sub target_region {
	my ($degradome_target, $genome, $output) = @_;
	system ("bowtie-build $genome $output\_bowtiedatabase");

	# create temporary file for bowtie alignments
	my $temp_bowtie = File::Temp->new();
	my $temp_rev_deg = File::Temp->new();
	system ("fastx_reverse_complement -i $degradome_target -o $temp_rev_deg");
	system ("bowtie -n 0 -l 10 -a $output\_bowtiedatabase -f $temp_rev_deg -S $temp_bowtie");
	# convert bowtie output to bedfile
	my $temp_bam = File::Temp->new();
	system ("samtools view -bS $temp_bowtie > $temp_bam");
	my $temp_bed = File::Temp->new();
	system ("bedtools bamtobed -i $temp_bam > $temp_bed");
	my $temp_alignments = File::Temp->new();
	open (BEDIN, $temp_bed) or die "can't open $!";
	open (BEDOUT, ">>$temp_alignments") or die "can't open $!";
	while (<BEDIN>) {
		if ($_ =~ /^(\S+\t)(\d+)\t(\d+)(\t\S+)(\t\d+)\t(.)/ ) {
			# make sure bed values are all positive
			if (($2 >= 20) and ($3 >= 20)) {
				if ($6 eq "-") { 
					my $value = $2-20;
					my $value2 = $2+20;
					print BEDOUT "$1$value\t$value2$4_$6_$2_$3$5\t$6\n";
					#print "$1$value\t$value2$4_$6_$2_$3$5\t$6\n";
				}
				elsif ($6 eq "+") {
					my $value = $3-20;
					my $value2 = $3+20;
					print BEDOUT "$1$value\t$value2$4_$6_$2_$3$5\t$6\n"
				}
			}
		}
	}
	
	# get fasta of output
	system ("bedtools getfasta -name -fi $genome -s -bed $temp_alignments -fo $output");
}

# take genome cleavage site and match sRNA to this region

sub cleavage_site_to_sRNA {
	my ($sRNA_targeting_degradome, $input_sequences, $mismatches, $output) = @_;
	#print "$output\n\n";
	my %degradome_hash;
	open (IN, $sRNA_targeting_degradome) or die "can't open $!";
	while (<IN>) {
		if ($_ =~ /([^-]+)-\d+\t(\w+)/) {
			exists($degradome_hash{$1}{$2}{$2});
		}
	}
	open (TARGET, $input_sequences) or die "can't open $!";
	open (OUT, ">>$output") or die "can't open $!";
	my $deg_seq = "";
	my %degradome_targets;
	my $degradome_hash_ref = \%degradome_hash; my $degradome_targets_ref = \%degradome_targets;
	while (<TARGET>) {
		if ($_ =~ /^>([^-]+)([^:]+::[^:]+)/) {
			if ($deg_seq eq $1) {
				my $coordinates = $2;
				my $target_seq = <TARGET>;
				exists ($degradome_targets{$deg_seq}{$target_seq}{"$deg_seq$coordinates"}{1});
			}
			else {
				my $new_deg_seq = $1; my $coordinates = $2;
				$degradome_targets_ref = \%degradome_targets;
				my $target_seq = <TARGET>;
				if (exists ($degradome_targets{$deg_seq})) {
					#					my $degradome_hash_ref = \%degradome_hash; my $degradome_targets_ref = \%degradome_targets;
					miRNA_to_target ($degradome_hash_ref, $degradome_targets_ref, $deg_seq, $mismatches, $output)
				}
				delete($degradome_targets{$deg_seq});
				$deg_seq = $new_deg_seq;
				exists ($degradome_targets{$deg_seq}{$target_seq}{"$deg_seq$coordinates"}{1});
			}
		}
	}
	#my $degradome_hash_ref = \%degradome_hash; my 
	$degradome_targets_ref = \%degradome_targets;
	miRNA_to_target ($degradome_hash_ref, $degradome_targets_ref, $deg_seq, $mismatches, $output);
}


sub miRNA_to_target {
	my ($degradome_hash_ref, $degradome_targets_ref, $deg_seq, $mismatches, $output) = @_;
	my $miRNA_temp = File::Temp->new();
	my %new_degradome = %{ $degradome_hash_ref };
	my %new_targets = %{ $degradome_targets_ref };
	open (OUTPUT, ">>$output") or die "can't open $!";
	foreach my $miRNA_seq (keys %{$new_degradome{$deg_seq}}) {
		print $miRNA_temp ">$miRNA_seq\n$miRNA_seq\n";
	}
	foreach my $target_sequence (keys %{$new_targets{$deg_seq}}) {
		my $target_temp = File::Temp->new();
		print $target_temp ">name\n$target_sequence";
		my $temp_patman_output = File::Temp->new();
		system ("patman -D $target_temp -P $miRNA_temp -e $mismatches -o $temp_patman_output");
		seek $temp_patman_output, 0, 0 or die "Seek $temp_patman_output failed: $!\n";
		while (<$temp_patman_output>) {
			$_ =~ /\S+\t(\S+)\t(\d+)\t(\d+)\t(.).*/;
			if ($4 eq "-") {
				if (20 <= $3 && 20 >= $2) {
					foreach my $header (keys %{$new_targets{$deg_seq}{$target_sequence}}) {
						print OUTPUT "$header\t$1\t$2\t$3\t$4\n";
					}
				}
			}
		}
	}
}



# take sRNA with putative targets and get their putative precursors

sub get_intersect {
	my ($hairpin_predictions, $degradome_targets, $maximum_targets, $output) = @_;
	open (DEG_OUT, $degradome_targets) or die "can't open $!";
	my %miRNA_predictions;
	my %degradome_predictions;
	open (OUT_DEG, ">>$output\_degradome") or die "can't open $!";
	open (OUT_MIR, ">>$output\_miRNA.fa") or die "can't open $!";
	while (<DEG_OUT>) {
		$_ =~ /^([^_]+)_(.)_(\d+)_(\d+)::(\w+)\t(\w+)\t\d+\t(\d+)/;
		my $cleavage = $7 - 20; my $cleavage2 = $cleavage + 1;
		$degradome_predictions{"$5\t$3\t$4\t$6\t0\t$2\t$1\n"}{$6} = "$cleavage\t$cleavage2";
		#		print OUT_DEG "$4\t$2\t$3\t$5\t0\t$1\n";
		$miRNA_predictions{$6}++;
	}
	print "\n\n$maximum_targets\n\n";
	open (IN_HAIRPIN, $hairpin_predictions) or die "can't open $!";
	open (OUT_HAIRPIN, ">>$output\_hairpin.bed") or die "can't open $!";
	while (<IN_HAIRPIN>) {
		if ($_ =~ /([^_]+)_\(([^)]+)\)_[^_]+_[^_]+_([^_]+)_([^_]+)_(.)/) {
		#		if ($_ =~ /([^_]+)_([^_]+)_[^_]+_[^_]+_([^_]+)_([^_]+)_(.)/) {
			if (exists ($miRNA_predictions{$1}) && ($miRNA_predictions{$1} <= $maximum_targets)) {
				print OUT_HAIRPIN "$2\t$3\t$4\t$1\t0\t$5\n";
			} 
		}
	}
	foreach my $deg (keys %degradome_predictions) {
		foreach my $miRNA (keys %{$degradome_predictions{$deg}}) {
			if ($miRNA_predictions{$miRNA} <= $maximum_targets) {
				my $degchomp = $deg;
				chomp $degchomp;
				print OUT_DEG "$degchomp\t$degradome_predictions{$deg}{$miRNA}\n";
			}
		}
	}
	foreach my $miRNA (keys %miRNA_predictions) {
		if ($miRNA_predictions{$miRNA} <= $maximum_targets) {
			print OUT_MIR ">$miRNA\n$miRNA\n";
		}
	}
}

1;
