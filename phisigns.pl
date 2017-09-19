#!/usr/bin/perl -w
#======================================================================================
## Author: Bhakti Dwivedi @ University of South Florida, FL, USA
## Description:  PhiSiGns Standalone version-1.0
##		         -Generates list of shared conserved genes ("Signature genes" or SiGs)
##                amongst the phage genomes as listed in the input file
##               -Generates FASTA file, and alignment for the selected SiG
##		         -Generates potential primer pairs for the selected SiG
## Usage: see README
#======================================================================================

use strict;
use Data::Dumper;

use Getopt::Long;
use File::Basename;
use List::Util qw(max);
#use Data::Dumper; 

use Bio::Seq;
use Bio::SeqIO;
use Bio::SeqFeature::Generic;
use Bio::DB::GenBank;
use Bio::DB::Fasta;
use Bio::Index::Fasta;
use Bio::Tools::Run::StandAloneBlast;
BEGIN { $ENV{CLUSTALDIR} = '/home/bin/clustalw/' }  # set CLUSTALW executable dir path
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Align::Utilities qw(aa_to_dna_aln);
use Bio::AlignIO;

#Environmental variables - expect them to be in the PATH
my $BLASTALL = 'blastall'; # set BLASTALL executable dir path
my $FORMATDB = 'formatdb'; # set FORMATDB executable dir path


#variables
my $EVALUE = 10; #blast e-value to search sequence similarity
my $COVERAGE = 10; #coverage (alignment length/avg len (query+hit) to search sequence similarity
my $MODE = 0;
my $DATA_DIR = 'data/';
my $OUTPUT_DIR = 'output/';
my $TMP = 'tmp/';
my $GENBANK = $TMP.'phage_gbk.txt';
my $BLASTOUT_EXT = '.blastoutput';
my $CDS_TABLE = $DATA_DIR.'cdstbl.ffd';
my $CDS_SEQS = $DATA_DIR.'cds_nucleotide.fa';
my $PPARAMS = 'pparams.txt';
my $IUPAC_DGN = { 'A' => 1,  # Adenine
                  'C' => 1,  # Cytosine
                  'G' => 1,  # Guanine
                  'T' => 1,  # Thymine
                  'M' => 2,  # A or C
                  'R' => 2,  # A or G
                  'W' => 2,  # A or T
                  'S' => 2,  # C or G
                  'Y' => 2,  # C or T
                  'K' => 2,  # G or T
                  'V' => 3,  # A or C or G
                  'H' => 3,  # A or C or T
                  'D' => 3,  # A or G or T
                  'B' => 3,  # C or G or T
                  'N' => 4   # G, A, T or C
		};
my $IUPAC_GC = {  'A' => 0,  # Adenine
                  'T' => 0,  # Thymine
				  'M' => 1,  # A or C
                  'R' => 1,  # A or G
                  'W' => 0,  # A or T
                  'S' => 2,  # C or G
                  'Y' => 1,  # C or T
                  'K' => 1,  # G or T
                  'V' => 2,  # A or C or G
                  'H' => 1,  # A or C or T
                  'D' => 1,  # A or G or T
                  'B' => 2,  # C or G or T
                  'N' => 2   # G, A, T or C
		};
my $IUPAC_COMPL = {'A' => 'T',
                  'C' => 'G',
                  'G' => 'C',
                  'T' => 'A',
		  		  'M' => 'K',
                  'R' => 'Y',
                  'W' => 'W',
                  'S' => 'S',
                  'Y' => 'R',
                  'K' => 'M',
                  'V' => 'B',
                  'H' => 'D',
                  'D' => 'H',
                  'B' => 'V',
                  'N' => 'N'
		};
our $WIN_SIZE = 10; #minium primer sequence length allowed from conserved region
our $STD_TMP = 37; #for the computation of Delta_G.
our $CONC_PRIMER = 200; #in nM change the concentration of primer - 2000000000 Convert from nanomoles to moles
our $CONC_SALT = 50; #in mM change the concentration of salt
our $CONC_MG = 0; #in mM change the concentration of salt
our $CONC_dNTP = 0; #in mM change the concentration of salt
our $TM_DIFF = 10; #max Tm differencs allowed
our $GC_DIFF = 10; #max gc difference allowed
our $N_NUM = 3; #number of 'N''s allowed in a primer sequence
our $GC_CLAMP = 0.60; #probability of 'G''s or 'C''s (not > 3) allowed within last 5 bases from the 3'end of primers.
our $MAX_3_STB = -9.0; #maximum 3' stability. DeltaG should not be more than -9kcal/mol for 5 bases from 3' end of primers
our $LOOP_BASES = 3; #number of min bases allowed in a loop when screening for hairpins,based on a ref
our $STEM_BASES = 4; #number of min bases allowed in a stem when screening for hairpins,based on a ref


#input
my $starttime = localtime;
my $command = basename($0)." ".join(' ',@ARGV)."\n";
my %params;
GetOptions(\%params,
	   'd=s', #database dir
	   'm=i', #mode
	   'aln_only', #display alignment only
	   'e=f', #evalue
	   'c=f', #coverage
	   'i=s', #input file
	   'g=s', #sig group number in the form: SiG_#
	   'a=s', #user uploaded alignment file
	   'id=s' #id of current job
	   );

#check if input is valid
if(exists $params{d}) {
	unless(-d $params{d}) {
		die "ERROR: database directory does not exists\n";
	}
} else {
    die "ERROR: you did not specify the database directory\n";
}
if(exists $params{m}) {
    unless($params{m} == 0 || $params{m} == 1) {
        die "ERROR: invalid value for mode. Mode 0 or 1 needs to be given! See README\n";
    }
} else {
    $params{m} = $MODE;
}
if(exists $params{e}) {
    unless($params{e} > 0) {
        die "ERROR: invalid value for evalue\n";
    }
} else {
    $params{e} = $EVALUE;
}
if(exists $params{c}) {
    unless($params{c} > 0) {
        die "ERROR: invalid value for coverage\n";
    }
} else {
    $params{c} = $COVERAGE;
}
if(exists $params{i}) {
	unless(-e $params{i}) {
		die "ERROR: cannot find input file\n";
	}
} elsif($params{m} == 0) {
    die "ERROR: you did not specify the input file\n";
}

#create sub directories
if (-d $OUTPUT_DIR) {
	&write_phisigns_log($params{id},"$OUTPUT_DIR already exists\n");
} else{
	mkdir ($OUTPUT_DIR,0777) || print $!;
}
if($params{m} == 0){
	if (-d $DATA_DIR) {
		&write_phisigns_log($params{id},"$DATA_DIR already exists\n");
	} else{
		mkdir ($DATA_DIR,0777) || print $!;
	}
	if (-d $TMP) {
		&write_phisigns_log($params{id},"$TMP already exists\n");
	} else{
		mkdir ($TMP,0777) || print $!;
	}
}


if($params{m} == 0) { #identify signature genes
	&write_phisigns_log($params{id},"Mode 0: Identifying signature genes in phages\n\n");
	&write_phisigns_log($params{id},"Start time Mode 0: $starttime\n");
	#read the user-selected phage genome accession #'s from the input file
	my %selected;
	open(IN, "<", $params{i}) or die "ERROR: could not open the input file: $! \n";
		while(<IN>) {
			chomp();
			$selected{$_} = 1;
		}
	close(IN);
	#check for number of phages
	die "STDERR Please select at least two phage genomes!\n\n" if(scalar(keys %selected) <= 1);

	#functions to get genbank flat files from NCBI, perform BLASTP search
	&get_genbank_data(\%selected);
	my ($gg,$geneid_len) = &parse_genbank_data();
	&format_db(\%selected);
	&precal_blast(\%selected);
	&merge_bl_files(\%selected,$geneid_len);

	#functions to identify signature gene groups
	my ($all_pairs,$evalue_op,$paralog_pairs,$evalue_pp) = &parse_blast($gg,\%selected,$params{d});
	my ($sig_grps) = &find_ortholog($all_pairs,$evalue_op);
	my ($sig_nops) = &collapse_grps($sig_grps);
	my $n = &compile_sig_genes($sig_nops, $gg, scalar(keys %selected));
	&write_phisigns_log($params{id},"End time Mode 0: ".(localtime)."\n");
	&write_phisigns_log($params{id},"\nMode 0: Identifying signature genes finished!\n");
	&write_phisigns_log($params{id},"Total signature genes (SiG's) identified = $n\n\n");

} elsif($params{m} == 1) { #design PCR primers pairs
	&write_phisigns_log($params{id},"Mode 1: Designing primer pairs for a selected signature gene\n\n");
	&write_phisigns_log($params{id},"Start time Mode 1: $starttime\n");
	my $align;
	if(exists $params{a}) {
		if($params{a} !~ m/.aln$/){
			die "ERROR: invalid clustalw file extension!\n";
			#exit(0);
		}
		$align = $params{a};
		&write_phisigns_log($params{id},"Using user uploaded alignment!\n");
	} else {
		if (-e $OUTPUT_DIR.$params{id}.'.sig_gps_detail.txt') {
			#get genes for selected SiG with program-generated alignment
			my $genes = &getSigGenes($OUTPUT_DIR.$params{id}.'.sig_gps_detail.txt',$params{g});
			&write_phisigns_log($params{id},"Using program generated alignment!\n");
			my $fasta = $OUTPUT_DIR.$params{id}.'.seqs';
			&createSeqsFile($genes,$fasta);
			#function to align selected SiG
			my ($prots,$seqs) = &dna_protein_alignment($fasta);
			$align = $OUTPUT_DIR.$params{id}.'.aln';
			&protein_alignment($align,$prots,$seqs);
		}
 		else {
			die "sig_gps_detail does not exist! Please execute phisigns in -m 0\n";
			#exit(0);
		}
	}
	#to display program generated alignment only
	if(exists $params{aln_only}) {
		&write_phisigns_log($params{id},"End time Mode 1: ".(localtime)."\n");
		&write_phisigns_log($params{id},"\nMode 1: Sequence alignment generated for the selected signature gene!\n\n");
		exit(0);
	} else {
		#Functions to design primers
		my $iupac_consensus = &make_consensus($align);
#		print STDERR Dumper $iupac_consensus;
		my $con_regions = &conserved_regions($iupac_consensus);
		my $pparams = &getPrimerParams($params{d}.$PPARAMS);
		&write_phisigns_log($params{id},"Searching for primer seqs in conserved regions\n");
		my ($region_windows,%pp,%uniq,%uniq_region_windows,%primers,%fwd,%rev,$primervals,$primervals_pdt);
		foreach my $key (keys %$con_regions) {
			$region_windows = &sliding_window($key,$con_regions->{$key});
		 		foreach my $window(keys %$region_windows){
					$uniq{$window}=$region_windows->{$window};
				}
			}
		&write_phisigns_log($params{id},"Searching for primer seqs in conserved regions finished\n");
		#search for primer pairs based on the parameter min and max values
		&write_phisigns_log($params{id},"Searching for primer pairs\n");
		$primervals = &primer_design(\%uniq,$pparams);
		$primervals_pdt = &get_pdt_length($primervals,\%pp,$pparams);
		&primer_design_2($primervals_pdt,\%primers,$pparams);
		my $pnum=0;
		if(scalar(keys %primers)) {
			my $pairs = $OUTPUT_DIR.$params{id}.'.pairs.txt';
			my $pairs_det = $OUTPUT_DIR.$params{id}.'.pairs_detail.txt';
			&get_unique_primers(\%primers,\%fwd,\%rev);
			$pnum = &print_primer_pairs(\%fwd,\%rev,$pairs,$pairs_det,$pparams);
		}
		&write_phisigns_log($params{id},"Searching for primer pairs finished\n");
   
	    &write_phisigns_log($params{id},"End time Mode 1: ".(localtime)."\n");
	    &write_phisigns_log($params{id},"\nMode 1: Designing primer pairs for selected signature genes finished!\n");
	    &write_phisigns_log($params{id},"Total primer pairs found = $pnum\n\n");
	}
}
#***************************************************************************************************************************
# Prints out the work progress to the log file
# Argument: comments
sub write_phisigns_log {
	my ($id,$comment) = @_;
	my $file = $OUTPUT_DIR.(defined $id ? $id : 'phisigns').'.log';
	open (LOG,">>$file") or die "couldn't open the $file file!";
	print LOG $comment if($comment);
	close(LOG);
}
#***************************************************************************************************************************
# Main PhiSiGns functions
#***************************************************************************************************************************
# Retreive genbank flat files from NCBI
# Arguments: array of selected refseq accession numbers
sub get_genbank_data {
	my $selected = shift;
	&write_phisigns_log($params{id},"Getting GenBank flat files for selected phages\n");
	my $writer = Bio::SeqIO->new( 
	    '-format' => 'genbank', 
	    '-file'   => ">$GENBANK" 
	); 
	foreach my $acc(keys %$selected){
		my $gb = Bio::DB::GenBank->new();
			 my $seq= $gb->get_Seq_by_acc($acc);
			    $writer->write_seq( $seq ); 
	}
	&write_phisigns_log($params{id},"Getting GenBank flat files for selected phages finished\n");
	return 0;
}
#***************************************************************************************************************************
# Parse GenBank files to retreive selected phage CDS info, CDS nucleotide FASTA, and others
# Arguments: None
sub parse_genbank_data{
	&write_phisigns_log($params{id},"Parsing GenBank flat files\n");
	my (%gg,%refseq,%geneid_len);
	my $reader = Bio::SeqIO->new( 
		   '-format' => 'genbank', 
		   '-file'   => "$GENBANK" 
	);
	# print the CDS table and CDS nucleotide file
	open (CDSTBL, ">", $CDS_TABLE) || die "couldn't open the $CDS_TABLE file!";
	open (CDS_NUCL, ">", $CDS_SEQS) || die "couldn't open the $CDS_SEQS file!";

	#parse genbank file, extract required info
	while (my $seq = $reader->next_seq ) {
		my $gacc =  $seq->accession_number;
	  	my @classification = $seq->species->classification;
		        pop @classification;
		        shift @classification;
			my $taxonomy = join(";", (@classification));
		my $phageclass;
		foreach my $classification (@classification){
			if ($classification =~ m/^\w+viridae$/i) {
				$phageclass = $classification;
			}
		}
		foreach my $ft ($seq->top_SeqFeatures) {
		  if ( $ft->primary_tag eq 'source' ) {  
			 my $cdscount=0;
			 my $glen = $ft->length;
			 my ($ggi) = $ft->has_tag('db_xref') && $ft->each_tag_value('db_xref');
			 my ($gmol) = $ft->has_tag('mol_type') && $ft->each_tag_value('mol_type');
			 my ($gorg) = $ft->has_tag('organism') && $ft->each_tag_value('organism');
			 	my @org_tag = split(/ /, $gorg);
			 	my $abr = substr($org_tag[0], 0, 3);
			 	my $gorg_tag = join("_", ($abr, $org_tag[$#org_tag]));
			 my ($phagehost) = $ft->has_tag('host') && $ft->each_tag_value('host');
			 	my @gid = split (/:/, $ggi);
			 	my $gid = $gid[1];
			      	my $writer_protein = Bio::SeqIO->new( 
			      		'-format' => 'fasta'); 
				if(exists $refseq{$gid}){
			      	$writer_protein = Bio::SeqIO->new( 
			      		'-format' => 'fasta', 
			      		'-file'   => ">>$TMP$gacc.fa"); 
				}
				else{
			      	$writer_protein = Bio::SeqIO->new( 
			      		'-format' => 'fasta', 
			      		'-file'   => ">$TMP$gacc.fa"); 
				}
				$refseq{$gid}=$gorg;
		     foreach my $f (grep { $_->primary_tag eq 'CDS' } $seq->get_SeqFeatures) {
				$cdscount++;
				my ($start,$end,$nclen,@start,@end);
			   	if ( $f->location->isa('Bio::Location::SplitLocationI'))  {
  				    my @sublocs = $f->location->sub_Location();
				    my $len=0;
				    #foreach my $location ( sort { $a->start <=> $b->start }  @sublocs ) {
				    foreach my $location ( $f->location->sub_Location) {
						push(@start,$location->start);
						push(@end,$location->end);
						$len+= abs(($location->start)-($location->end));
						$len+=1;
				    }
					$start = join(":",@start);
					$end = join(":",@end);
					$nclen = $len;
			  	} else{
					$start=$f->start;
					$end=$f->end;
					$nclen = $f->length;
				}
				my ($gname);
					if ( $f->has_tag('gene') ) {($gname) = $f->each_tag_value('gene');
					} 
					elsif ($f->has_tag('product') ) {($gname) = $f->each_tag_value('product');
					}
       				my $aalen = int(($nclen/3)-1);
				my $strand = $f->strand;
				my ($locus) = $f->has_tag('locus_tag') && $f->each_tag_value('locus_tag');
				my ($ref) = $f->has_tag('protein_id') && $f->each_tag_value('protein_id');
				my ($gi)  = ($f->has_tag('db_xref') && $f->each_tag_value('db_xref')) || ($f->has_tag('protein_id') && $f->each_tag_value('protein_id'));

				my ($proid)  = $f->has_tag('product') && $f->each_tag_value('product');
				#deprekated
				#my @gi = split (/:/, $gi);
				#$gi = $gi[1];
				$gg{$ref}=$gacc;
				$geneid_len{$ref}=$aalen;
		            	#push (@{$gg{$gacc}},$ref);
				my ($translation) = $f->has_tag('translation') && $f->each_tag_value('translation');
		       		unless( $gi && $ref && $gname && $translation ) {
			   		&write_phisigns_log($params{id},"STDERR not fully annotated CDS $gi,$ref,$gname,skipping\n");
					 next;
				}
				#deprekated
				#my $cds_obj = $f->spliced_seq($translation);
				my $cds_obj = $f->spliced_seq();
			     	my $dna = $cds_obj->seq;
				print CDS_NUCL "$ref:$dna\n";
				my $proteinseq = Bio::PrimarySeq->new(
						'-seq' => $translation,
					     	'-display_id' => sprintf("$ref")
				);
				$writer_protein->write_seq( $proteinseq );
				print CDSTBL join ("\t", $ref, $gi, $nclen, $aalen, $start, $end,$strand, $locus, $proid, $gorg, $phageclass)."\n";
		     	}
		   }
		 }     
	      }
	close (CDSTBL);
	close (CDS_NUCL);
	&write_phisigns_log($params{id},"Parsing GenBank flat files finished\n");
	return (\%gg,\%geneid_len);
}
# Execute FORMATDB
# Arguments: selected phages refseq accession number
sub format_db {
	my $selected = shift;
	&write_phisigns_log($params{id},"Formatting database using FORMATDB\n");
	foreach my $key(keys %$selected){
		system ("$FORMATDB -i $TMP$key.fa -p T -o T -n $TMP$key");
	}
	&write_phisigns_log($params{id},"Formatting database using FORMATDB finished\n");
	return 0;
}
# Execute BLASTALL
# Arguments: selected phages refseq accession number
sub precal_blast {
	&write_phisigns_log($params{id},"BLASTP all against all search\n");
	my $selected = shift;
	foreach my $key_1 (keys %$selected){
		foreach my $key_2(keys %$selected){
			my $file = $TMP.$key_1.".".$key_2.$BLASTOUT_EXT;
			if($key_1 eq $key_2){
		          system ("$BLASTALL -p blastp -i $TMP$key_1.fa -d $TMP$key_2 -e $EVALUE -o $file -m 8 -v 2 -b 2");}
			if($key_1 ne $key_2){
		          system ("$BLASTALL -p blastp -i $TMP$key_1.fa -d $TMP$key_2 -e $EVALUE -o $file -m 8 -v 1 -b 1");}
		}
       }
	&write_phisigns_log($params{id},"BLASTP all against all search finished\n");
	return 0;
}
# Parse BLASTP comparison files into BLASTP output per phage
# Argument: blast file names generated in precal_blast function
sub merge_bl_files{
	my ($selected,$geneid_len) = @_;
	&write_phisigns_log($params{id},"Parsing blast files per phage genome\n");
	foreach my $key_1(keys %$selected){
	open (MERGE, ">", $DATA_DIR.$key_1.$BLASTOUT_EXT)|| die "couldn't open $DATA_DIR.$key_1.$BLASTOUT_EXT file!";
		foreach my $key_2(keys %$selected){
		     my $file = $TMP.$key_1.".".$key_2.$BLASTOUT_EXT;
		     open (DATA, "<", $file)|| die "couldn't open $file file!";
			 while (my $line = <DATA>){
				chomp $line;
				my @cols = split(/\t/, $line);
					my $avg = ($geneid_len->{$cols[0]} + $geneid_len->{$cols[1]})/2;
					my $coverage = ($cols[3]/$avg)*100;
		     	 	print MERGE join ("\t", $cols[0],$cols[1],sprintf("%.2f",$coverage),$cols[10])."\n";
		   	 }
	      	     close (DATA);
		}
	close (MERGE);
	}
	&write_phisigns_log($params{id},"Parsing blast files per phage genome finished\n");
	clean_bl_files();
	return 0;
}
# Parse blast results for user-selected phage genomes
# Arguments: array of gene ids, selected phage genome, directory
sub parse_blast {
	my ($gg,$selected,$dir) = @_;
	&write_phisigns_log($params{id},"Parsing blast results\n");
	my (@args,$paralog,%paralog_ids,$parsed,$all);
	$parsed = 0;
	$all = scalar(keys %$selected)+1;
	my (%ortholog_pairs,%evalue,%paralog_pairs,%evalue_paralogpairs,$file,$tmp);
	foreach my $gid_query (keys %$selected) {
		$file = $dir.$DATA_DIR.$gid_query.$BLASTOUT_EXT;
		if(-e $file) {
			open(BLAST, "<", $file) or die "ERROR: couldn't open $file file: $! \n";
			while(<BLAST>) {
				chomp();
				@args = split(/\t/);
				if(exists $gg->{$args[1]} && $params{e} >= $args[3] && $params{c} <= $args[2] && $args[0] ne $args[1]) {
					$tmp = $args[0].' '.$args[1];
					if($gid_query eq $gg->{$args[1]}) { #paralog: genome id query == genome id subject
						unless(exists $paralog_ids{$args[0]}){
							$paralog_pairs{(sort($args[0],$args[1]))[0]}->{(sort($args[0],$args[1]))[1]}++;
							$evalue_paralogpairs{$tmp} = $args[3].' '.$args[2];
							$paralog_ids{$args[0]} = 1;
						}
					} elsif(!exists $evalue{$tmp}) {
						$ortholog_pairs{(sort($args[0],$args[1]))[0]}->{(sort($args[0],$args[1]))[1]}++;
						$evalue{$tmp} = $args[3].' '.$args[2];
					}
				}
			}
			close (BLAST);
		} else {
			die "ERROR: can not find BLAST output file for $gid_query in $dir$DATA_DIR\n";
		}
		$parsed++;
	}

	&write_phisigns_log($params{id},"Parsing blast results finished\n");
	return (\%ortholog_pairs,\%evalue,\%paralog_pairs,\%evalue_paralogpairs);
}
# This is for the ortholog identification
# Arguments: $ortholog_pairs and their E-values from parse_blast function
sub find_ortholog {
	my ($ortholog_pairs,$evalue) = @_;
	&write_phisigns_log($params{id},"Searching for shared genes\n");
	my (%hash,@sig_grps,@implied,@second,@tmp,$check,$count);
	foreach my $first (sort keys %{$ortholog_pairs}) {
		$count = $check = 0;
		@tmp = @implied = @second = ();
		foreach my $second (keys %{$ortholog_pairs->{$first}}) {
			next if($ortholog_pairs->{$first}->{$second} < 2);
			unless($check){
				push(@tmp,$first);
				$check = 1;
			}
			push(@tmp,$second);
			$count++;
			push(@implied,$second);
                }
		@implied = sort @implied;
		foreach my $i (0..$#implied) {
			foreach my $j ($i+1..$#implied) {
				if (exists $ortholog_pairs->{$implied[$i]}->{$implied[$j]}) {
					delete $ortholog_pairs->{$implied[$i]}->{$implied[$j]};
					$count++;}
			}
		}
		next if !$count;
		push(@sig_grps,[@tmp]);
	}
	write_phisigns_log($params{id},"Searching for shared genes finished\n");
	return (\@sig_grps);
}
# This is to collapse signature gene groups that share a gene.
# Arguments: the sig_grps array from the ortholog identification function
sub collapse_grps {
	my $sig_grps_ref = shift;
	&write_phisigns_log($params{id},"Combining related signature gene groups\n");
	my (%gs_grp,@sig_nets,$ids,%gids,$net,$check);
	foreach my $i (0..scalar(@$sig_grps_ref)-1) {
		foreach my $gid (@{$sig_grps_ref->[$i]}) {
			push(@{$gs_grp{$gid}}, $i);
		}
	}
	foreach my $i (0..scalar(@$sig_grps_ref)-1) {
		next unless($sig_grps_ref->[$i]);
		$ids = $sig_grps_ref->[$i];
		delete $sig_grps_ref->[$i];
		%gids = ();
		foreach my $gid (@$ids) {
			next if(exists $gids{$gid});
			$gids{$gid} = 1;
			foreach my $j (@{$gs_grp{$gid}}) {
				next unless($sig_grps_ref->[$j]);
				$net = $sig_grps_ref->[$j];
				delete $sig_grps_ref->[$j];
				push(@$ids, @$net);
			}
		}
		push(@sig_nets,join(' ', sort keys %gids));
	}
	&write_phisigns_log($params{id},"Combining related signature gene groups finished\n");
	return (\@sig_nets);
}
# Put together the signature gene list
# Arguments: $sig_nops from collapse function, geneids, and # of selected phage genomes
sub compile_sig_genes {
	my ($sig_nops, $gg, $num) = @_;
	&write_phisigns_log($params{id},"Putting together signature genes\n");
	my ($sigs,$sigs_detail,%pp_count,@present,$presentn,$present_percentage,%ids,@tmp,$count,%counts,$mean,%funcs);
	$count = 0;
	foreach my $cols (@$sig_nops) {
		@tmp = split(/\s+/, $cols);
		foreach(@tmp) {
			$ids{$_} = 1;
			$pp_count{$count}{$gg->{$_}}++;
			$counts{$count}->{$_} = 1;
		}
		$count++;
	}
	&get_siginfo_cds(\%ids);

	$count = 1;
	$sigs = $sigs_detail = '';
	foreach my $c (sort {scalar(keys %{$pp_count{$b}}) <=> scalar(keys %{$pp_count{$a}})} keys %pp_count) {
		@tmp = keys %{$counts{$c}};
		$presentn = scalar(keys %{$pp_count{$c}});
		push(@present, $presentn);
		$present_percentage = ($num > 0 ? sprintf("%.2f",$presentn*100/$num) : 0);
		$mean = 0;
		%funcs = ();
		foreach my $id (@tmp) {
#			print STDERR Dumper $id,$ids{$id},@tmp,$count,$presentn,$present_percentage unless(@{$ids{$id}} && $count && $presentn && $present_percentage);
			$sigs_detail .= join("\t", 'SiG_'.$count,$presentn,$present_percentage,@{$ids{$id}})."\n";
			$mean += $ids{$id}->[1];
			$funcs{$ids{$id}->[0]}++;
		}
		$sigs_detail .= "\n";
		$sigs .= join("\t", 'SiG_'.$count,$presentn,$present_percentage,(sort {$funcs{$b} <=> $funcs{$a}} keys %funcs)[0],int($mean/scalar(@tmp)))."\n"; #get the mean amino acid length per signature group
		$count++;
	}
	if($count == 1) {
		&write_phisigns_log($params{id},"\nSTDERR No signature gene identified in the phages selected!\n");
		&write_phisigns_log($params{id},"STDERR Either relax the E-value, Coverage cut-off or change the selection of phages!\n\n");
	} else {
		&write_phisigns_log($params{id},"Putting together signature genes finished\n");
	}
	&print_sig_genes($sigs, $sigs_detail);
	return scalar(@$sig_nops);
}
# Print the signature gene list
# Arguments: @sorted_sigs and @sorted_sig_detail from sort_sig_genes function
sub print_sig_genes {
	my ($sorted_sigs, $sorted_sigs_detail) = @_;
	&write_phisigns_log($params{id},"Printing the list of signature genes\n");

	open(SIG_GPS, ">",$OUTPUT_DIR.$params{id}.'.sig_gps.txt') or die "ERROR: could not open the sig_gps file: $! \n";
	print SIG_GPS $sorted_sigs;
	close (SIG_GPS);
	open(SIG_DETAIL, ">",$OUTPUT_DIR.$params{id}.'.sig_gps_detail.txt') or die "ERROR: could not open the sig_gps_detail file: $! \n";
	print SIG_DETAIL $sorted_sigs_detail;
	close (SIG_DETAIL);

	&write_phisigns_log($params{id},"Printing the list of signature genes finished\n");
}
# Converts the dna sequences into amino acid
# Arguments: selected SiG sequence fasta file
sub dna_protein_alignment {
	my $fasta = shift;
	&write_phisigns_log($params{id},"Translating dna -> protein for alignment\n");
	my $seqin = Bio::SeqIO->new('-format' => 'fasta', '-file' => $fasta);
	my (%seqs,@prots,$pseq);
	while(my $seq = $seqin->next_seq) {
		$seqs{$seq->display_id} = $seq;
		my $protein = $seq->translate(-codontable_id => 11); #bacterial codon table
		$pseq = $protein->seq();
		$pseq =~ s/\*$//;
		if($pseq =~ /\*/) {
			&write_phisigns_log($params{id},"\tprovided a cDNA (".$seq->display_id.") sequence with a stop codon\n");
			#exit(0);
			$pseq =~ s/\*//g;
		}
		$protein->seq($pseq);
		push(@prots,$protein);
	}
	if(scalar(@prots) < 2 ) {
		&write_phisigns_log($params{id},"\tNeed at least 2 cDNA sequences to proceed\n");
		exit(0);
	}
	&write_phisigns_log($params{id},"Translating dna -> protein for alignment finished\n");
	return (\@prots,\%seqs);
}
# CLUSTALW on protein sequences and generate alignments
# Arguments: alignment file, aa and nt sequences for selected SiG
sub protein_alignment {
	my ($align,$prots,$seqs) = @_;
	&write_phisigns_log($params{id},"Aligning protein sequences\n");
	my @params = ('-ktuple' => 2, '-matrix' => 'BLOSUM'); #parameters can be changed
	my $aln_factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
	my $aa_aln = $aln_factory->align($prots);
	my $dna_aln = aa_to_dna_aln($aa_aln, $seqs);
	my @each = $dna_aln->each_seq();
	$aa_aln->map_chars('\.','-'); #clustal writes '.' instead of '-'
	#my $new_dna_aln = $dna_aln->remove_gaps('-'); # remove by position
	my $alnout = Bio::AlignIO->new('-format' => 'clustalw', '-file'   => ">$align");
	#$alnout->write_aln($aa_aln);#to print the protein alignment
	$alnout->write_aln($dna_aln);#to print the nucleotide alignment
	&write_phisigns_log($params{id},"Aligning protein sequences finished\n");
	return 0;
}
# Create the consensus IUPAC string from multiple nucleotide sequence alignment
# Arguments: alignment file
sub make_consensus {
	my $align = shift;
	&write_phisigns_log($params{id},"Making IUPAC consensus for the selected $align\n");
        my $aio = Bio::AlignIO->new('-format' => 'clustalw', '-file' => $align);
	my $dna_aln = $aio->next_aln;
	my $iupac_consensus = $dna_aln->consensus_iupac();
	$iupac_consensus =~ tr/a-z/./; #"." denotes gaps in the sequence alignment
	&write_phisigns_log($params{id},"Making IUPAC consensus for the selected SiG finished\n");
	return $iupac_consensus;
}
# Retrieve the conserved blocks from the consensus iupac string (>=100bp)
# Arguments: IUPAC consensus string
sub conserved_regions {
	my $consensus = shift;
	&write_phisigns_log($params{id},"Extracting conserved regions from consensus IUPAC\n");
	if(length($consensus) < 100){
		&write_phisigns_log($params{id},"STDERR The consensus sequence length should be more than 100bp!\n");
		return 0;
	}
	my %con_regions;
	LOOP: for my $str ( 0 .. length( $consensus ) - 20 ) {
		my( $match ) = $consensus =~ m[.{$str}([ACGT].{0,18}[ACGT])];
		if(!defined $match){next LOOP;}
		m[$match] and next LOOP for keys %con_regions;
		my $start = index($consensus, $match);
		my $end = ($start+length $match)-1;
		my ($add_start,$add_end,$len);
			#incase when the conserved residue falls at the start and end of the string
			if($start ==0){	$add_start = $start;
					$add_end = $end+8;
					$len = ($add_end-$add_start)+1;
			}elsif($start ==1){$add_start = $start-1;
					   $add_end = $end+7;
					   $len = ($add_end-$add_start)+1;
			}elsif($start ==2){$add_start = $start-2;
					   $add_end = $end+6;
					   $len = ($add_end-$add_start)+1;
			}elsif($start ==3){$add_start = $start-3;
					   $add_end = $end+5;
					   $len = ($add_end-$add_start)+1;
			}elsif($end ==(length ($consensus)-1)){$add_start = $start-8;
							       $add_end = $end;
							       $len = ($add_end-$add_start)+1;
			}elsif($end == length ($consensus)-2){$add_start = $start-7;
							      $add_end = $end+1;
							      $len = ($add_end-$add_start)+1;
			}elsif($end == length ($consensus)-3){$add_start = $start-6;
							      $add_end = $end+2;
							      $len = ($add_end-$add_start)+1;
			}elsif($end == length ($consensus)-4){$add_start = $start-5;
							      $add_end = $end+3;
							      $len = ($add_end-$add_start)+1;
			}else{$add_start = $start-4;
			      $add_end = $end+4;
			      $len = ($add_end-$add_start)+1;
			     }
		my $cons_reg = substr($consensus, $add_start, $len);
		my $sp = $add_start+1;
		my $ep = $add_end+1;
		my $position = "$sp".":"."$ep";
		if($cons_reg =~ m/\.+/){# exclude regions with gaps
			next;
		}
		$con_regions{$cons_reg}=$position; #Check: same conserved region will overwrite older one
	}
	&write_phisigns_log($params{id},"Extracting conserved regions from consensus IUPAC finished\n");
	if(scalar keys (%con_regions) == 1) {
		&write_phisigns_log($params{id},"STDERR Only one conserved region identified!\n");
		&write_phisigns_log($params{id},"STDERR Unable to design primer pairs!\n");
		exit(0);
	}
	&write_phisigns_log($params{id},"Extracting conserved regions from consensus IUPAC finished\n");
	return \%con_regions;
}
# Get primer parameters minimum and maximum values
# Arguments: primer parameter file 'pparams.txt'
sub getPrimerParams {
	my ($file) = @_;
	my %pparams;
	my $data = &getIdFileData($file);
	foreach my $s
('minl','maxl','mingc','maxgc','minbtm','maxbtm','minstm','maxstm','minnnh','maxnnh','mindg','minpl','maxpl','dgn','gcclamp','max3stb','compl')
{
	    $pparams{$s} = $data->{$s}||'';
	}
	return \%pparams;
}
# sliding window for each conserved region identified
# Arguments: conserved_regions and position from conserved_regions function
sub sliding_window {
	my ($con_regions, $position) = @_;
	my (%windows,$window,$sp,$ep);
	my @tmp = split(/\:/,$position);
	my $start = $tmp[0];
	my $end = $tmp[1];
	for(my $winsize = length($con_regions); $winsize >= $WIN_SIZE; $winsize--) {
		for(my $i = 0; $i <= length($con_regions)-$winsize; $i++) {
			$sp = $start;
			$ep = $end;
	  		$window = substr($con_regions,$i,$i+$winsize);
			$sp += $i;
			$ep = $sp + length($window) - 1;
			$windows{$window} = $sp.':'.$ep;
		}
	}
	return \%windows;
}
# Obtain primer sequences that qualify the basic primer parameter filters
# Arguments: primer sequences and parameters minima and maxima
sub primer_design {
	my ($region_windows,$pparams) = @_;
	my @primervals=();
	my ($start,$stop,$condition,$count,$len,$GC,$degeneracy,$basicTm,$SATm,$primer_min,$primer_max);
	my ($GC_min,$GC_max,$basicTm_min,$basicTm_max,$SATm_min,$SATm_max);
	foreach my $window (keys %$region_windows) {
		my @position = split(/\:/, $region_windows->{$window});
		$start = $position[0];
		$stop = $position[1];
		$condition = &condition_check($window);
		if($condition eq 'deoxy') {
			$len = length $window;
			next if($len < $pparams->{minl} || $len > $pparams->{maxl});
			$GC = &GC($condition,$window);
			next if($GC < $pparams->{mingc} || $GC > $pparams->{maxgc});
			$degeneracy = 1;
			$basicTm = &basic_Tm($window);
			next if($basicTm < $pparams->{minbtm} || $basicTm > $pparams->{maxbtm});
			$SATm = &salt_adjusted_Tm($window);
			next if($SATm < $pparams->{minstm} || $SATm > $pparams->{maxstm});
			push(@primervals,[$window,$condition,$len,$start,$stop,sprintf("%.2f",$GC),sprintf("%.2f",$GC),$degeneracy,sprintf("%.2f",$basicTm),sprintf("%.2f",$basicTm),sprintf("%.2f",$SATm),sprintf("%.2f",$SATm)]);
		} elsif($condition eq 'iupac') {
			next if($window =~ m/^N/ || $window =~ m/N$/);
			$count = ($window =~ tr/N//);
			next if($count > $N_NUM);
			$primer_min = &primer_min($window);
			$primer_max = &primer_max($window);
#			print STDERR Dumper "primer",$primer_min,$primer_max;
			$len = length $window;
#			print STDERR Dumper "len",$len;
			next if($len < $pparams->{minl} || $len > $pparams->{maxl});
			$GC_min = &GC($condition, $window);
			$GC_max = &GC($condition, $window);
#			print STDERR Dumper "GC",$GC_min,$GC_max;
			next if($GC_min < $pparams->{mingc} || $GC_max > $pparams->{maxgc});
			$degeneracy = &degeneracy($window);
#			print STDERR Dumper "degen",$degeneracy,$pparams->{dgn};
			next if($degeneracy > $pparams->{dgn});
			$basicTm_min = &basic_Tm($primer_min);
			$basicTm_max = &basic_Tm($primer_max);
#			print STDERR Dumper "basicTm",$basicTm_min,$basicTm_max;
			next if($basicTm_min < $pparams->{minbtm} || $basicTm_max > $pparams->{maxbtm});
			$SATm_min = &salt_adjusted_Tm($primer_min);
			$SATm_max = &salt_adjusted_Tm($primer_max);
#			print STDERR Dumper "SATm",$SATm_min,$SATm_max;
			next if($SATm_min < $pparams->{minstm} || $SATm_max > $pparams->{maxstm});
			push(@primervals,[$window,$condition,$len,$start,$stop,sprintf("%.2f",$GC_min),sprintf("%.2f",$GC_max),$degeneracy,sprintf("%.2f",$basicTm_min),sprintf("%.2f",$basicTm_max),sprintf("%.2f",$SATm_min),sprintf("%.2f",$SATm_max)]);
		}
	}
	return (\@primervals);
}
# check the product length-to determine the forward and reverse primers
# Arguments: primer sequence array filtered through the primer_design function, parameters minima and maxima
sub get_pdt_length{
	my ($primervals,$pp,$pparams) = @_;
	my ($pdt,$fdp,$rvp,$values,%primervals_pdt);
	foreach my $p1 (@$primervals) {
		for my $p2 (@$primervals) {
		        #check the product length
			$pdt = abs($p1->[3] - $p2->[4]) + 1;  # start of one - end of the other
			if($pdt >= $pparams->{minpl} && $pdt <= $pparams->{maxpl}) { #product length filter
				$fdp = $rvp = '';
				#determine the forward and reverse primers
				if($p1->[3] < $p2->[3]) {
					$fdp = $p1->[0];
					$rvp = &rev_complement($p2->[0]);
					$pp->{$fdp}->{$rvp}=1;
					$values = join("\t", 'FOR',$fdp,$p1->[1],$p1->[2],$p1->[3],$p1->[4],$p1->[5],$p1->[6],$p1->[7],$p1->[8],$p1->[9],$p1->[10],$p1->[11]);
					$primervals_pdt{$fdp}=$values;
					#push(@primervals_pdt,['FOR',$fdp,$p1->[1],$p1->[2],$p1->[3],$p1->[4],$p1->[5],$p1->[6],$p1->[7],$p1->[8],$p1->[9],$p1->[10],$p1->[11]]);
					$values = join("\t", 'REV',$rvp,$p2->[1],$p2->[2],$p2->[3],$p2->[4],$p2->[5],$p2->[6],$p2->[7],$p2->[8],$p2->[9],$p2->[10],$p2->[11]);
					$primervals_pdt{$rvp}=$values;
					#push(@primervals_pdt,['REV',$rvp,$p2->[1],$p2->[2],$p2->[3],$p2->[4],$p2->[5],$p2->[6],$p2->[7],$p2->[8],$p2->[9],$p2->[10],$p2->[11]]);
				} else {
					$fdp = $p2->[0];
					$rvp = &rev_complement($p1->[0]);
					$pp->{$fdp}->{$rvp}=1;
					$values = join("\t", 'FOR',$fdp,$p2->[1],$p2->[2],$p2->[3],$p2->[4],$p2->[5],$p2->[6],$p2->[7],$p2->[8],$p2->[9],$p2->[10],$p2->[11]);
					$primervals_pdt{$fdp}=$values;
					#push(@primervals_pdt,['FOR',$fdp,$p2->[1],$p2->[2],$p2->[3],$p2->[4],$p2->[5],$p2->[6],$p2->[7],$p2->[8],$p2->[9],$p2->[10],$p2->[11]]);
					$values = join("\t", 'REV',$rvp,$p1->[1],$p1->[2],$p1->[3],$p1->[4],$p1->[5],$p1->[6],$p1->[7],$p1->[8],$p1->[9],$p1->[10],$p1->[11]);
					$primervals_pdt{$rvp}=$values;
					#push(@primervals_pdt,['REV',$rvp,$p1->[1],$p1->[2],$p1->[3],$p1->[4],$p1->[5],$p1->[6],$p1->[7],$p1->[8],$p1->[9],$p1->[10],$p1->[11]]);
				}
			}
		}
	}
	return (\%primervals_pdt);
}
# Filter primer pairs obtained from get_pdt_length function
# Arguments: primer pairs from get_pdt_length function, parameters minima and maxima
sub primer_design_2 {
	my ($primervals_pdt,$primers,$pparams) = @_;
	my ($molwt,%nucsCount,$H,$S,$M_dG,$NNTm,$gcclamp,$dG3,$count);
	my ($primer_seq,$primer_min,$primer_max,$len,$molwt_min,$molwt_max,@IUpairVals_min,@IUpairVals_max,$base0,$base,$MdG_min,$MdG_max,$H_min,$H_max,$S_min,$S_max,$NNTm_min,$NNTm_max,$gcclamp_min, $gcclamp_max,$dG3_min,$dG3_max,@temp,$H_endsmin,$S_endsmin,$H_endsmax,$S_endsmax);
	  foreach my $key(keys %{$primervals_pdt}) {
		    my @primer = split(/\t/, $primervals_pdt->{$key});
			  if($primer[2] eq 'deoxy'){
				%nucsCount = ();
				foreach my $nucs (qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)) {
					$nucsCount{$nucs} = &countNeighbors($nucs,$primer[1]);
				}
				$H = &deltaH(\%nucsCount, $primer[1]);
				$S = &deltaS(\%nucsCount, $primer[1]);
				$M_dG = &recal_deltaG($H,$S);
				next if($M_dG < $pparams->{mindg});
				$NNTm = &nearestNeighborTM($H, $S);
				next if($NNTm < $pparams->{minnnh} || $NNTm > $pparams->{maxnnh});
				if($pparams->{gcclamp} eq 'y') {
					$gcclamp = &GC_clamp($primer[1]);
					next if($gcclamp > $GC_CLAMP);
				}
				if($pparams->{max3stb} eq 'y') {
					$dG3 = &max_3_stability($primer[1]);
					next if($dG3 < $MAX_3_STB);
				}
				my $values = join("\t",$primer[0],$primer[3],$primer[4],$primer[5],$primer[6],$primer[7],$primer[8],$primer[9],$primer[10],$primer[11],$primer[12],sprintf("%.2f",$NNTm),sprintf("%.2f",$NNTm),sprintf("%.2f",$M_dG),sprintf("%.2f",$M_dG));
				$primers->{$primer[1]} = $values;
			} elsif($primer[2] eq 'iupac') {
				@IUpairVals_min = @IUpairVals_max = (0.0,0.0,0.0); #deltaG, deltaH, deltaS
				$len = length($primer[1]);
				for(my $i=0; $i<$len; $i++) {
              				$base0 = substr($primer[1],$i,1);
              				$base = substr($primer[1],$i+1,1);
					@temp = &calcIUpair($base0, $base, $i, 'min',$primer[1]);
					for(my $j=0; $j<3; $j++) {
						$IUpairVals_min[$j]+=$temp[$j];
					}
					@temp = &calcIUpair($base0, $base, $i, 'max',$primer[1]);
					for(my $j=0; $j<3; $j++) {
						$IUpairVals_max[$j]+=$temp[$j];
					}
				}
				$H_min = $IUpairVals_min[1];
				$S_min = $IUpairVals_min[2];
				($H_endsmin,$S_endsmin) = &HS_correct_ends($primer[1],'min');
				$H_min += $H_endsmin;
				$S_min += $S_endsmin;
				$H_max = $IUpairVals_max[1];
				$S_max = $IUpairVals_max[2];
				($H_endsmax,$S_endsmax) = &HS_correct_ends($primer[1],'max');
				$H_max += $H_endsmax;
				$S_max += $S_endsmax;
				my $salt_S = &deltaS_correction($primer[1]);
				$S_min += $salt_S;
				$S_max += $salt_S;
				$MdG_min = &recal_deltaG($H_min,$S_min);
				$MdG_max = &recal_deltaG($H_max,$S_max);
				next if($MdG_min < $pparams->{mindg});
				$NNTm_max = &nearestNeighborTM($H_min,$S_min);
				$NNTm_min = &nearestNeighborTM($H_max,$S_max);
				next if($NNTm_min < $pparams->{minnnh} || $NNTm_max > $pparams->{maxnnh});
#				print STDERR Dumper "NNTm",$NNTm_min,$NNTm_max;
				if($pparams->{gcclamp} eq 'y') {
					$gcclamp = &GC_clamp($primer[1]);
					next if($gcclamp > $GC_CLAMP);
				}
				if($pparams->{max3stb} eq 'y') {
					$dG3 = &max_3_stability($primer[1]);
					next if($dG3 < $MAX_3_STB);
				}
#				print STDERR Dumper "found";
				my $values = join("\t",$primer[0],$primer[3],$primer[4],$primer[5],$primer[6],$primer[7],$primer[8],$primer[9],$primer[10],$primer[11],$primer[12],sprintf("%.2f",$NNTm_min),sprintf("%.2f",$NNTm_max),sprintf("%.2f",$MdG_min),sprintf("%.2f",$MdG_max));
				$primers->{$primer[1]} = $values;
			}
		}
}
#Prints unique forward and reverse primer pairs
#Arguments: primer pairs and their info from primer_design subroutines
sub get_unique_primers {
	my ($primers,$fwd,$rev) = @_;
	&write_phisigns_log($params{id},"Define forward & reverse primers\n");
	my ($fwd_count, $rev_count)=(0,0);
	foreach my $p (keys %$primers){
		my @p_val = split(/\t/, $primers->{$p});
		if($p_val[0] eq 'FOR'){
			$fwd_count++;
			my $self_dimer = &primer_dimer("Fwd_$fwd_count",$p,$p);
			my $hairpin = &primer_hairpin("Fwd_$fwd_count",$p);
			my $values = join("\t",$p,$p_val[1],$p_val[2],$p_val[4],$p_val[5],$p_val[6],$p_val[7],$p_val[8],$p_val[9],$p_val[10],$p_val[11],$p_val[12],$p_val[13],$p_val[14],$self_dimer,$hairpin);
			$fwd->{$fwd_count} = $values;
		}elsif($p_val[0] eq 'REV'){
			$rev_count++;
			my $self_dimer = &primer_dimer("Rev_$rev_count",$p,$p);
			my $hairpin = &primer_hairpin("Rev_$rev_count",$p);
			my $values = join("\t",$p,$p_val[1],$p_val[2],$p_val[4],$p_val[5],$p_val[6],$p_val[7],$p_val[8],$p_val[9],$p_val[10],$p_val[11],$p_val[12],$p_val[13],$p_val[14],$self_dimer,$hairpin);
			$rev->{$rev_count} = $values;
		}
	}
	if($fwd_count == 0 || $rev_count == 0) {
		&write_phisigns_log($params{id},"\nSTDERR Not enough Forward or Reverse to make a primer pair!\n");
		&write_phisigns_log($params{id},"STDERR Relax your parameters!\n\n");
		#exit(0);
	}
	&write_phisigns_log($params{id},"Defining forward & reverse finished\n");
	return 0;
}
#Prints forward and reverse primer pairs
#Arguments: selected forward and reverse.
sub print_primer_pairs {
	my ($fwd,$rev,$pairs,$pairs_det,$pparams) = @_;
	&write_phisigns_log($params{id},"Printing primer pairs\n");
	my ($pdt,$gc_diff,$tm_diff_BT,$tm_diff_SAT,$tm_diff_NNT,$tm_mismatch,$cross_dimer,$count,@p2);
	$count = 0;
	#open the output files
	open(PAIRS,">$pairs") or die "ERROR: could not open file: $! \n";
	open(PAIRS_DETAIL,">$pairs_det") or die "ERROR: could not open file: $! \n";
	#header row
	print PAIRS '#'.join("\t", "Pair","5'-Forward-3'","Len","Degen","5'-Reverse-3'","Len","Degen","Pdt Len","Tm Mismatch","Complementarity")."\n";
	print PAIRS_DETAIL '#'.join("\t", "Primer", "5'-seq -3'","Len","Start","Degen","GC% Min/Max","BasicTm Min/MAx","SaltAdjTm Min/Max","NNTm Min/Max","Delta G Min/Max","Pdt Len","Tm Mismatch","Self-dimer","Cross-dimer","hairpin")."\n";

	foreach my $fwd_p (keys %$fwd){
		my @p1=split(/\t/, $fwd->{$fwd_p});
		foreach my $rev_p (keys %$rev){
			my $compl_check;
			@p2 = split(/\t/, $rev->{$rev_p});
			$pdt = abs($p1[2] - ($p2[2]+$p2[1])+1);   # start of one - (start of the other+length of the other)
			$gc_diff = abs($p1[3] - $p2[3]);   # min GC of one and GC Tm of other
			$tm_diff_BT = abs($p1[6] - $p2[6]);   # min Tm of one and min Tm of other
			$tm_diff_SAT = abs($p1[8] - $p2[8]);  # min Tm of one and min Tm of other
			$tm_diff_NNT = abs($p1[10] - $p2[10]);  # min Tm of one and min Tm of other

			$tm_mismatch = ($tm_diff_BT+$tm_diff_SAT+$tm_diff_NNT)/3;  #take the avg tm diff over three Tm

			if($tm_diff_BT <= $TM_DIFF && $tm_diff_SAT <= $TM_DIFF && $tm_diff_NNT <= $TM_DIFF && $gc_diff <= $GC_DIFF && $pdt >= $pparams->{minpl} && $pdt <= $pparams->{maxpl}) {
				$cross_dimer = &primer_dimer("$fwd_p.$rev_p",$p1[0],$p2[0]);
				#output only if Pass cross-dimer check
				if($pparams->{compl} eq 'y'){
					if($cross_dimer eq 'Pass' && $p1[14] eq 'Pass' && $p2[14] eq 'Pass' && $p1[15] eq 'Pass' && $p2[15] eq 'Pass') {
						$compl_check = "Pass";
						$count++;
						print PAIRS join("\t",$count,$p1[0],$p1[1],$p1[5],$p2[0],$p2[1],$p2[5],$pdt,sprintf("%.2f",$tm_mismatch),$compl_check)."\n";

						print PAIRS_DETAIL join ("\t", "SIGPR_".$count."_FOR",$p1[0],$p1[1],$p1[2],$p1[5],$p1[3].'/'.$p1[4],$p1[6].'/'.$p1[7],$p1[8].'/'.$p1[9],$p1[10].'/'.$p1[11],$p1[12].'/'.$p1[13], $pdt,sprintf("%.2f",$tm_mismatch),$p1[14],$cross_dimer,$p1[15])."\n";
						print PAIRS_DETAIL join ("\t", "SIGPR_".$count."_REV",$p2[0],$p2[1],$p2[2],$p2[5],$p2[3].'/'.$p2[4],$p2[6].'/'.$p2[7],$p2[8].'/'.$p2[9],$p2[10].'/'.$p2[11],$p2[12].'/'.$p2[13], $pdt,sprintf("%.2f",$tm_mismatch),$p2[14],$cross_dimer,$p2[15])."\n";
						print PAIRS_DETAIL "\n";
						}
					}
				#output all possibilies if compl unchecked
				if($pparams->{compl} eq 'n'){
					if($cross_dimer eq 'Fail' && $p1[14] eq 'Fail' && $p2[14] eq 'Fail' && $p1[15] && 'Fail' && $p2[15] eq 'Fail') {
							$compl_check = "Fail";
						}
					elsif($cross_dimer eq 'Pass' && $p1[14] eq 'Pass' && $p2[14] eq 'Pass' && $p1[15] && 'Pass' && $p2[15] eq 'Pass') {
							$compl_check = "Pass";
						}
					else{
							$compl_check = "Warning";
						}
						$count++;
						print PAIRS join("\t",$count,$p1[0],$p1[1],$p1[5],$p2[0],$p2[1],$p2[5],$pdt,sprintf("%.2f",$tm_mismatch),$compl_check)."\n";

						print PAIRS_DETAIL join ("\t", "SIGPR_".$count."_FOR",$p1[0],$p1[1],$p1[2],$p1[5],$p1[3].'/'.$p1[4],$p1[6].'/'.$p1[7],$p1[8].'/'.$p1[9],$p1[10].'/'.$p1[11],$p1[12].'/'.$p1[13], $pdt,sprintf("%.2f",$tm_mismatch),$p1[14],$cross_dimer,$p1[15])."\n";
						print PAIRS_DETAIL join ("\t", "SIGPR_".$count."_REV",$p2[0],$p2[1],$p2[2],$p2[5],$p2[3].'/'.$p2[4],$p2[6].'/'.$p2[7],$p2[8].'/'.$p2[9],$p2[10].'/'.$p2[11],$p2[12].'/'.$p2[13], $pdt,sprintf("%.2f",$tm_mismatch),$p2[14],$cross_dimer,$p2[15])."\n";
						print PAIRS_DETAIL "\n";
					}
			}
		}
	}
	close(PAIRS);
	close(PAIRS_DETAIL);
	if($count == 0) {
		unlink($pairs);
		unlink($pairs_det);
		&write_phisigns_log($params{id},"\nSTDERR No Primer Pairs found! This could be because\n");
		&write_phisigns_log($params{id},"\tSTDERR primers have high complementarity within itself or within primer pairs\n");
		&write_phisigns_log($params{id},"\tSTDERR primers show presence of secondary structure such as hairpins\n");
		&write_phisigns_log($params{id},"\tSTDERR primer parameters are too stringent\n");
		&write_phisigns_log($params{id},"\tSTDERR selected sequences are too divergent\n");
		#exit(0);
	}
	&write_phisigns_log($params{id},"Printing primer pairs finished\n");
	return $count;
}
#***************************************************************************************************************************
# Primer parameters computation functions
#***************************************************************************************************************************
#calculate the GC content (in percentage) in the given sequence
#Arguments: sequence string
sub GC {
	my ($condition, $sequence) = @_;
	my ($count_GC,$GC_content) = 0;
	if ($condition eq 'deoxy') {
		$count_GC = &count_GC($sequence);
	}elsif($condition eq 'iupac') {
		my @residues = split(//,$sequence);
		foreach my $residue (@residues){
	 		if($residue =~ m/G|C/){
	   			$count_GC++;
	    		} else{
	 			$count_GC += $IUPAC_GC->{$residue} / $IUPAC_DGN->{$residue};
			}
     		}
	}
	$GC_content = ($count_GC / length($sequence) * 100);
	return $GC_content;
}
# to compute the degeneracy of the sequence
# Arguments: sequence string
sub degeneracy {
        my $sequence=shift;
        my $dgn = 1;
        for (my $i = 0; $i < length($sequence); $i++) {
               my $base = substr($sequence, $i, 1);
               $dgn *= $IUPAC_DGN->{$base};
        }
	return $dgn;
}
#Basic Tm Calculations.
#Arguments: sequence string
sub basic_Tm {
	my $sequence = shift;
	my $seq_len = length($sequence);
	my $n_AT = &count_AT($sequence);
	my $n_CG = &count_GC($sequence);
	if ($seq_len < 14) {
	        return (2 * ($n_AT) + 4 * ($n_CG));
	} else {
	        return (64.9 + 41*(($n_CG-16.4)/$seq_len));
	}
}
#Salt-Adjusted Tm Calculations
#Arguments: sequence string
sub salt_adjusted_Tm {
	my $sequence = shift;
	my $primer_len = length($sequence);
	my $n_AT = count_AT($sequence);
	my $n_CG = count_GC($sequence);
	my $r = &log10($CONC_SALT/1000);
	if($primer_len < 14) {
		return (2 * ($n_AT) + 4*($n_CG) - 16.6*&log10(0.05) + 16.6*$r);
	} else {
		#return (81.5  +  41*($n_CG/$primer_len) - (500/$primer_len) + 16.6*$r - 0.62*F;
		#return (81.5  +  16.6*($r)  +  41*($n_CG/$primer_len) - (675/$primer_len));
		return (100.5 + 41.0*($n_CG/$primer_len)-(820.0/$primer_len)+16.6*$r);
	}
}
#GC_CLAMP calculation.  Not more than 3 G's or C's (or 0.60 prob) within the last five bases from the 3' end of primers
#Arguments: last 5 bases from 3' end of forward and reverse primers obtained from get_3Pri_5bases function
sub GC_clamp{
	my $sequence = shift;
	my $last5bases = &get_3Pri_5bases ($sequence);
	my @last5bases = split(//,$last5bases);
	my $prob=0;
	foreach my $i (0 .. scalar @last5bases-1){
		if($last5bases[$i] eq 'G' || $last5bases[$i] eq 'C'){$prob+=1;}
		elsif(exists $IUPAC_GC->{$last5bases[$i]}){
			$prob+=($IUPAC_GC->{$last5bases[$i]}/$IUPAC_DGN->{$last5bases[$i]});
		}
	}
	return sprintf("%.2f",$prob/5);
}
#Computation of maximum 3' stability.  DeltaG kcal/mol using equation: dG = dH -TdS should not be less than -9kcal/mol
#Arguments: last 5 bases from 3' end of forward and reverse primers
sub max_3_stability {
	my $sequence = shift;
	my $last5bases = &get_3Pri_5bases ($sequence);
	my ($H,$S,$dG,%nucsCount);
	my ($base0,$base,$H_min,$H_max,$S_min,$S_max,$salt_S,$dG_min,$dG_max,@temp,@IUpairVals_min,@IUpairVals_max);
	my $condition = &condition_check($last5bases);
	if($condition eq 'deoxy') {
		foreach my $nucs (qw(AA AC AG AT CA CC CG CT GA GC GG GT TA TC TG TT)) {
			$nucsCount{$nucs} = &countNeighbors($nucs,$last5bases);
		}
		$H = &deltaH(\%nucsCount, $last5bases);
		$S = &deltaS(\%nucsCount, $last5bases);
		$dG = &recal_deltaG($H,$S);
	} else {
		for(my $i=0; $i<length $last5bases; $i++) {
			$base0 = substr($last5bases,$i,1);
			$base = substr($last5bases,$i+1,1);
			@temp = &calcIUpair($base0, $base, $i, 'min',$last5bases);
			for(my $j=0; $j<3; $j++) {
				$IUpairVals_min[$j]+=$temp[$j];
			}
			@temp = &calcIUpair($base0, $base, $i, 'max',$last5bases);
			for(my $j=0; $j<3; $j++) {
				$IUpairVals_max[$j]+=$temp[$j];
			}
		}
		$H_min = $IUpairVals_min[1];
		$H_max = $IUpairVals_max[1];
		my ($H_ends,$S_ends) = &HS_correct_ends ($last5bases,'min');
		$H_min+=$H_ends;
		$S_min+=$S_ends;
		$H_max = $IUpairVals_max[1];
		$S_max = $IUpairVals_max[2];
		($H_ends,$S_ends) = &HS_correct_ends ($last5bases,'max');
		$H_max+=$H_ends;
		$S_max+=$S_ends;
		$salt_S = &deltaS_correction($last5bases);
		$S_min += $salt_S;
		$S_max += $salt_S;
		$dG_min = &recal_deltaG($H_min,$S_min);
		$dG_max = &recal_deltaG($H_max,$S_max);

		$dG = ($dG_min < $dG_max ? $dG_min : $dG_max);
	}
	return $dG;
}
#Computation of DeltaG kcal/mol using equation: dG = dH -TdS; used in the current delta G computation.
#Arguments: $deltaH, $deltaS, and std temperature in degree celcius
sub recal_deltaG {
	my ($Delta_H, $Delta_S) = @_;
	my $Delta_G = 0;
	$Delta_G = $Delta_H-((273.15 + $STD_TMP)*($Delta_S/1000)); #dG = dH-TdS
	return $Delta_G;
}
#Computation of Nearest Neighbor Tm
#Arguments: deltaH and deltaS
sub nearestNeighborTM {
	my ($Delta_H, $Delta_S) = @_;
	my $NN_TM = 0.00;
	$NN_TM = (($Delta_H * 1000) / ($Delta_S + (1.987 * log($CONC_PRIMER/4000000000)))) - 273.15;
	return $NN_TM;
}
#check the primers for dimer formation, sliding algorithm-aligns each possible combination
#Arguments: primer; fwd vs fwd; rvp: rvp or fwd vs rvp
sub primer_dimer {#self-dimer and cross-dimer
	my ($count,$fp,$rp_53) = @_;
	my $fpl = length($fp);
	my $rpl = length($rp_53);
	my $rp = reverse($rp_53);
	my($lg_match,$max_match,$max_index,$min_pd_deltaG) = (0,0,0,0);
	my ($new_binding,$max_binding,$string,$status,$dimer);
	$status = 0;

	#slide the forward primer on the reverse primer
	#stable 3'end extensible primer dimer
	for(my $i=$fpl-1; $i>0; $i--) {
		my $first = substr($fp,$i,$fpl-$i);
		my $second = substr($rp,0,$rpl-$i);
		my ($match,$binding) = &compare($first,$second);
		($lg_match,$string,$new_binding) = &find_longest_string($first,$binding);
		next unless $lg_match > 0;#if no matched bases or less than 3 consecutive bases, do not report
		#my $pd_deltaG = &deltaG_primer_dimer($string);
		#if($pd_deltaG < $min_pd_deltaG){#store the worse-case deltaG - to draw hairpin
		#search for the dimer with maximum number of match bases
		if($lg_match > $max_match) {
			#$min_pd_deltaG=$pd_deltaG;
			$max_match = $lg_match;
#			$max_binding = $new_binding;
#			$max_index = $i;
			$status = 1;
			last;
		}
	}
	if($status) {
		$dimer = 'Warning';
		#print DIM "---------------------------------\n";
		#print DIM "3' extensible primer-dimers \n";
		#print DIM "---------------------------------\n";
#			print DIM '#'.$count, "\n";
#			print DIM "Best Matches: ", $max_match, "\n";
#			print DIM "Delta G: ", $min_pd_deltaG, " kcal/mole","\n\n";
#			my @pattern = split(//,$max_binding);
#			my $spaces = " " x $max_index;
#			print DIM "\t5' ",$fp, "\n";
#			print DIM "\t   ",$spaces;
#				foreach my $j(@pattern){
#					if($j=='2'){
#						print DIM "|";
#					}elsif($j=='1'){
#						print DIM ":";
#					}else{
#						print DIM " ";
#					}
#				}
#			print DIM "\n";
#			print DIM "\t3' ",$spaces,$rp,"\n\n";
		return $dimer;
	} else {
		$dimer = 'Pass';
	}
	($lg_match,$max_match,$max_index,$min_pd_deltaG) = (0,0,0,0);#reinitialize these variables
	$max_binding = '';
	$status = 0;

	#slide the reverse primer on the forward primer
	#stable non-extensible primer dimer; simply reduces primer conc in solution
	for(my $j=$rpl-1; $j>=0; $j--) {
		my $first = substr($rp,$j,$rpl-$j);
		my $second = substr($fp,0,$fpl-$j);
		#check for complementarity
		my ($match,$binding)=&compare ($first,$second);
		($lg_match,$string,$new_binding)=&find_longest_string($first,$binding);
		next unless $lg_match > 0;
		#my $pd_deltaG = &deltaG_primer_dimer($string);
		#if($pd_deltaG < $min_pd_deltaG){#store the worse-case deltaG - to draw hairpin
		#search for the dimer with maximum number of match bases
		if($lg_match > $max_match) {
			#$min_pd_deltaG=$pd_deltaG;
			$max_match=$lg_match;
#			$max_binding=$new_binding;
#			$max_index=$j;
			$status = 1;
			last;
		}
	}
	if($status) {
		$dimer = 'Warning'; 
		#print DIM "---------------------------------\n";
		#print DIM "Non extensible primer-dimers \n";
		#print DIM "---------------------------------\n";
#			print DIM '#'.$count, "\n";
#			print DIM "Best Matches: ", $max_match, "\n";
#			print DIM "Delta G: ", $min_pd_deltaG, " kcal/mole","\n";
#			my @pattern = split(//,$max_binding);
#			my $spaces = " " x $max_index;
#			print DIM "\t5' ",$spaces, $fp, "\n";
#			print DIM "\t   ",$spaces;
#				foreach my $j(@pattern){
#					if($j=='2'){
#						print DIM "|";
#					}elsif($j=='1'){
#						print DIM ":";
#					}else{
#						print DIM " ";
#					}
#				}
#			print DIM "\n";
#			print DIM "\t3' ",$rp,"\n\n";
		return $dimer;
	} else {
		$dimer = 'Pass';
	}
	#print DIM '#'.$count, "\n";
	#print DIM "None Found!\n\n";
	return $dimer;
}
#check for secondary structure- hairpin formation
#Arguments:primer sequence (fwd or rev)
sub primer_hairpin {
	my ($count,$primer) = @_;
	my $pl=length($primer);
	my ($lg_match,$max_match,$hp_deltaG,$min_hp_deltaG) = (0,0,0,0);
	my (@pattern,$new_binding,$binding,$max_binding,$string,$stem1,$stem2,$max_stem1,$max_stem2,$loop,$status,$hairpin);
	$status = 0;
	foreach my $i (0..($pl - 1 - $LOOP_BASES - 2 * $STEM_BASES)) { #loop through the primer sequence
		$stem1 = reverse(substr($primer, 0, $STEM_BASES+$i)); #stem length should be 4 bases and above
		$stem2 = substr($primer, $STEM_BASES+$i+$LOOP_BASES);
		my ($match,$binding) = &compare($stem1,$stem2);
		($lg_match,$string,$new_binding) = &find_longest_string($primer,$binding);
#		next unless($lg_match > 0);#if no matched bases or < 3 consecutive bases, do not report
		#$hp_deltaG = &deltaG_primer_dimer($string);
		#if($hp_deltaG < $min_hp_deltaG){#store the worse-case deltaG - to draw hairpin
		#search for the dimer with maximum number of match bases
		if($lg_match > $max_match) {
#			$min_hp_deltaG=$hp_deltaG;
			$max_match = $lg_match;
#			$max_binding=$new_binding;
#			$max_stem1 = $stem1;
#			$max_stem2 = $stem2;
#			$loop = substr($primer, $STEM_BASES+$i,$LOOP_BASES);
			$status = 1;
			last;
		}
	}

	if($status) {
		$hairpin = 'Warning'; #is there a way you can have this in red?
		#Following output opens a new tab, when Warning is clicked!
		#displays worse case scenario only based on # of matches, there could be more.
		#not sure to compute or show delta G for that matter.
		#print HP "---------------------------------\n";
		#print HP "Hairpin Formation\n";
		#print HP "---------------------------------\n";
#			print HP $primer,"\n";
#			print HP  '#'.$count, "\n";
#			print HP "Best Matches: ", $max_match, "\n";
#			print HP "Delta G: ", $min_hp_deltaG, " kcal/mole","\n";
#			print HP "  ",substr($loop,0,1)," \\\n";
#			print HP " /   ",$max_stem1,"\n";
#			print HP substr($loop,1,1),"    ";
#			@pattern = split(//,$max_binding);
#				foreach my $j(@pattern){
#					if($j=='2'){
#						print HP "|";
#					}elsif($j=='1'){
#						print HP ":";
#					}else{
#						print HP " ";
#					}
#				}
#			print  HP "\n";
#			print  HP " \\   ",$max_stem2,"\n";
#			print  HP "  ",substr($loop,2,1)," /\n\n";
	} else {
		$hairpin = 'Pass';
		#print HP '#'.$count, "\n";
		#print HP "None Found!\n\n";
	}
	return $hairpin;
}
#count the number of occurences of a character pair, used when string is dna
#Arguments: neighboring nt pairs, and sequence string
sub countNeighbors{
	my ($string,$sequence) = @_;
	my $count = 0;
	my $result = index($sequence, $string);
	while($result != -1) {
		$count++;
		$result = index($sequence, $string, $result+1);
	}
	return $count;
}
#The following computation are done using the Santa Lucia PNAS 1998 parameter values
#Computation of DeltaH kcal/mol, used when sequence string is deoxy
#Arguments: neighbor nt pairs, sequence string
sub deltaH {
	my ($pairVals,$sequence) = @_;
	my $Delta_H = 0.0;
		$Delta_H+= -7.9*$pairVals->{AA};
		$Delta_H+= -7.9*$pairVals->{TT};
		$Delta_H+= -7.2*$pairVals->{AT};
		$Delta_H+= -7.2*$pairVals->{TA};
		$Delta_H+= -8.5*$pairVals->{CA};
		$Delta_H+= -8.4*$pairVals->{GT};
		$Delta_H+= -7.8*$pairVals->{CT};
		$Delta_H+= -8.2*$pairVals->{GA};
		$Delta_H+=-10.6*$pairVals->{CG};
		$Delta_H+= -9.8*$pairVals->{GC};
		$Delta_H+= -8.0*$pairVals->{GG};
		$Delta_H+= -8.0*$pairVals->{CC};
		$Delta_H+= -8.5*$pairVals->{TG};
		$Delta_H+= -8.2*$pairVals->{TC};
		$Delta_H+= -8.4*$pairVals->{AC};
		$Delta_H+= -7.8*$pairVals->{AG};

		my $firstnucleotide = substr($sequence,0,1);
		if($firstnucleotide eq 'G' || $firstnucleotide eq 'C') {
			$Delta_H += 0.1;
		} elsif($firstnucleotide eq 'A' || $firstnucleotide eq 'T') {
			$Delta_H += 2.3;
		}
		my $lastnucleotide = substr($sequence,length($sequence)-1,1);
		if($lastnucleotide eq 'G' || $lastnucleotide eq 'C') {
			$Delta_H += 0.1;
		} elsif($lastnucleotide eq 'A' || $lastnucleotide eq 'T') {
			$Delta_H += 2.3;
		}
	return $Delta_H;
}
#Computation of DeltaS cal/k.mol,used when sequence string is deoxy
#Arguments: neighbor nt pairs, sequence string
sub deltaS {
	my($pairVals,$sequence) = @_;
	my $Delta_S = 0.0; 
		$Delta_S+= -22.2*$pairVals->{AA};
		$Delta_S+= -22.2*$pairVals->{TT};
		$Delta_S+= -20.4*$pairVals->{AT};
		$Delta_S+= -21.3*$pairVals->{TA};
		$Delta_S+= -22.7*$pairVals->{CA};
		$Delta_S+= -22.4*$pairVals->{GT};
		$Delta_S+= -21.0*$pairVals->{CT};
		$Delta_S+= -22.2*$pairVals->{GA};
		$Delta_S+= -27.2*$pairVals->{CG};
		$Delta_S+= -24.4*$pairVals->{GC};
		$Delta_S+= -19.9*$pairVals->{GG};
		$Delta_S+= -19.9*$pairVals->{CC};
		$Delta_S+= -22.7*$pairVals->{TG};
		$Delta_S+= -22.2*$pairVals->{TC};
		$Delta_S+= -22.4*$pairVals->{AC};
		$Delta_S+= -21.0*$pairVals->{AG};

		my $firstnucleotide = substr($sequence,0,1);
		if($firstnucleotide eq 'G' || $firstnucleotide eq 'C') {
			$Delta_S += -2.8;
		} elsif($firstnucleotide eq 'A' || $firstnucleotide eq 'T') {
			$Delta_S += 4.1;
		}
		my $lastnucleotide = substr($sequence,length($sequence)-1,1);
		if($lastnucleotide eq 'G' || $lastnucleotide eq 'C') {
			$Delta_S += -2.8;
		} elsif($lastnucleotide eq 'A' || $lastnucleotide eq 'T') {
			$Delta_S += 4.1;
		}
	my $salt_dS = &deltaS_correction ($sequence);
	$Delta_S += $salt_dS;
	return $Delta_S;
}
#Salt correction for deltaS
#Arguments: sequence string
sub deltaS_correction{
	my $sequence = shift;
	my $dS;
	my $salt_correction=$CONC_SALT/1000;
	$dS += 0.368 * (length($sequence) - 1) * log($salt_correction);
	return $dS;
}
#Computation of DeltaG kcal/mol.  Not used in the currect dG computation
#Arguments: neighbor nt pairs, sequence string
sub deltaG {
	my ($pairVals, $sequence) = @_;
	my $Delta_G = 0.0;
		$Delta_G+= -1.00*$pairVals->{AA};
		$Delta_G+= -1.00*$pairVals->{TT};
		$Delta_G+= -0.88*$pairVals->{AT};
		$Delta_G+= -0.58*$pairVals->{TA};
		$Delta_G+= -1.45*$pairVals->{CA};
		$Delta_G+= -1.44*$pairVals->{GT};
		$Delta_G+= -1.28*$pairVals->{CT};
		$Delta_G+= -1.30*$pairVals->{GA};
		$Delta_G+= -2.17*$pairVals->{CG};
		$Delta_G+= -2.24*$pairVals->{GC};
		$Delta_G+= -1.84*$pairVals->{GG};
		$Delta_G+= -1.84*$pairVals->{CC};
		$Delta_G+= -1.45*$pairVals->{TG};
		$Delta_G+= -1.30*$pairVals->{TC};
		$Delta_G+= -1.44*$pairVals->{AC};
		$Delta_G+= -1.28*$pairVals->{AG};

		my $firstnucleotide = substr($sequence,0,1);
		if($firstnucleotide eq 'G' || $firstnucleotide eq 'C') {
			$Delta_G += 0.98;
		} elsif($firstnucleotide eq 'A' || $firstnucleotide eq 'T') {
			$Delta_G += 1.03;
		}
		my $lastnucleotide = substr($sequence,length($sequence)-1,1);
		if($lastnucleotide eq 'G' || $lastnucleotide eq 'C') {
			$Delta_G += 0.98;
		} elsif($lastnucleotide eq 'A' || $lastnucleotide eq 'T') {
			$Delta_G += 1.03;
		}
	return $Delta_G;
}
#calculating neighboring iupac pairs, used when string is iupac
#Arguments: Base0, base, index, choice, and sequence string
sub calcIUpair {
	my ($base0, $base, $i, $choice,$seq) = @_;
	my (@IUpacBase0,@IUpacBase);
	my @reValue = (0,0,0);

	if(&is_base($base0) && &is_base($base)) { #for a condition as AA
		my $tmp = &getPairVal($base0.$base);
		@reValue = @$tmp;
	} elsif(&is_iupac_base($base0) && &is_base($base)) {#for a condition BA
		if($base0 eq 'M') {
			@IUpacBase0 = ('A','C');
		} elsif($base0 eq 'R') {
			@IUpacBase0 = ('A','G');
		} elsif($base0 eq 'W') {
			@IUpacBase0 = ('A','T');
		} elsif($base0 eq 'S') {
			@IUpacBase0 = ('C','G');
		} elsif($base0 eq 'Y') {
			@IUpacBase0 = ('C','T');
		} elsif($base0 eq 'K') {
			@IUpacBase0 = ('G','T');
		} elsif($base0 eq 'V') {
			@IUpacBase0 = ('A','C','G');
		} elsif($base0 eq 'H') {
			@IUpacBase0 = ('A','C','T');
		} elsif($base0 eq 'D') {
			@IUpacBase0 = ('A','G','T');
		} elsif($base0 eq 'B') {
			@IUpacBase0 = ('C','G','T');
		} elsif($base0 eq 'N') {
			@IUpacBase0 = ('A','C','G','T');
		}
		my $j = 0;
		foreach my $base0 (@IUpacBase0) {
			my $tmp = &getPairVal($base0.$base);
			if($j==0){
				@reValue = @$tmp;
			} else {
				foreach my $k (0..2) {
					if($choice eq 'min' && $reValue[$k] > $tmp->[$k]) {
						$reValue[$k] = $tmp->[$k];
					} elsif($choice eq 'max' && $reValue[$k] < $tmp->[$k]) {
						$reValue[$k] = $tmp->[$k];
					}
				}
			} $j++;
		}
	} elsif(&is_base($base0) && &is_iupac_base($base)) {#for a condition AB
		if($base eq 'M') {
			@IUpacBase = ('A','C');
		} elsif($base eq 'R') {
			@IUpacBase = ('A','G');
		} elsif($base eq 'W') {
			@IUpacBase = ('A','T');
		} elsif($base eq 'S') {
			@IUpacBase = ('C','G');
		} elsif($base eq 'Y') {
			@IUpacBase = ('C','T');
		} elsif($base eq 'K') {
			@IUpacBase = ('G','T');
		} elsif($base eq 'V') {
			@IUpacBase = ('A','C','G');
		} elsif($base eq 'H') {
			@IUpacBase = ('A','C','T');
		} elsif($base eq 'D') {
			@IUpacBase = ('A','G','T');
		} elsif($base eq 'B') {
			@IUpacBase = ('C','G','T');
		} elsif($base eq 'N') {
			@IUpacBase = ('A','C','G','T');
		}
		my $j = 0;
		foreach my $base (@IUpacBase) {
			my $tmp = &getPairVal($base0.$base);
			if($j==0){
				@reValue = @$tmp;
			} else {
				foreach my $k (0..2) {
					if($choice eq 'min' && $reValue[$k] > $tmp->[$k]) {
						$reValue[$k] = $tmp->[$k];
					} elsif($choice eq 'max' && $reValue[$k] < $tmp->[$k]) {
						$reValue[$k] = $tmp->[$k];
					}
				}
			} $j++;
		}
	} elsif(&is_iupac_base($base0) && &is_iupac_base($base)) {	
			if($base0 eq 'M') {
				@IUpacBase0 = ('A','C');
			} elsif($base0 eq 'R') {
				@IUpacBase0 = ('A','G');
			} elsif($base0 eq 'W') {
				@IUpacBase0 = ('A','T');
			} elsif($base0 eq 'S') {
				@IUpacBase0 = ('C','G');
			} elsif($base0 eq 'Y') {
				@IUpacBase0 = ('C','T');
			} elsif($base0 eq 'K') {
				@IUpacBase0 = ('G','T');
			} elsif($base0 eq 'V') {
				@IUpacBase0 = ('A','C','G');
			} elsif($base0 eq 'H') {
				@IUpacBase0 = ('A','C','T');
			} elsif($base0 eq 'D') {
				@IUpacBase0 = ('A','G','T');
			} elsif($base0 eq 'B') {
				@IUpacBase0 = ('C','G','T');
			} elsif($base0 eq 'N') {
				@IUpacBase0 = ('A','C','G','T');
			}
			
			if($base eq 'M') {
				@IUpacBase = ('A','C');
			} elsif($base eq 'R') {
				@IUpacBase = ('A','G');
			} elsif($base eq 'W') {
				@IUpacBase = ('A','T');
			} elsif($base eq 'S') {
				@IUpacBase = ('C','G');
			} elsif($base eq 'Y') {
				@IUpacBase = ('C','T');
			} elsif($base eq 'K') {
				@IUpacBase = ('G','T');
			} elsif($base eq 'V') {
				@IUpacBase = ('A','C','G');
			} elsif($base eq 'H') {
				@IUpacBase = ('A','C','T');
			} elsif($base eq 'D') {
				@IUpacBase = ('A','G','T');
			} elsif($base eq 'B') {
				@IUpacBase = ('C','G','T');
			} elsif($base eq 'N') {
				@IUpacBase = ('A','C','G','T');
			}
		my $j = 0;
		foreach my $base0 (@IUpacBase0) {
			foreach my $base (@IUpacBase) {
				my $tmp = &getPairVal($base0.$base);
				if($j==0){
					@reValue = @$tmp;
				} else {
					foreach my $k (0..2) {
						if($choice eq 'min' && $reValue[$k] > $tmp->[$k]) {
							$reValue[$k] = $tmp->[$k];
						} elsif($choice eq 'max' && $reValue[$k] < $tmp->[$k]) {
							$reValue[$k] = $tmp->[$k];
						}
					}
				} $j++;
			}
		}
	}
	return @reValue;
}
sub getPairVal { #get the nearest-neighbor values
	my $pair = shift;
	my @tmp;
	if($pair eq 'AA' || $pair eq 'TT') {
		@tmp = (-1.0,-7.9,-22.2);
	} elsif($pair eq 'AT') {
		@tmp = (-0.88,-7.2,-20.4);
	} elsif($pair eq 'TA') {
		@tmp = (-0.58,-7.2,-21.3);
	} elsif($pair eq 'CA') {
		@tmp = (-1.45,-8.5,-22.7);
	} elsif($pair eq 'TG') {
		@tmp = (-1.44,-8.5,-22.7);
	} elsif($pair eq 'GT') {
		@tmp = (-1.44,-8.4,-22.4);
	} elsif($pair eq 'AC') {
		@tmp = (-1.45,-8.4,-22.4);
	} elsif($pair eq 'CT') {
		@tmp = (-1.28,-7.8,-21.0);
	} elsif($pair eq 'AG') {
		@tmp = (-1.30,-7.8,-21.0);
	} elsif($pair eq 'GA') {
		@tmp = (-1.30,-8.2,-22.2);
	} elsif($pair eq 'TC') {
		@tmp = (-1.28,-8.2,-22.2);
	} elsif($pair eq 'CG') {
		@tmp = (-2.17,-10.6,-27.2);
	} elsif($pair eq 'GC') {
		@tmp = (-2.24,-9.8,-24.4);
	} elsif($pair eq 'GG' || $pair eq 'CC') {
		@tmp = (-1.84,-8.0,-19.9);
	}
	return \@tmp;
}
#terminal corrections for iupac strings using probability
#Arguments: sequence string, choice
sub HS_correct_ends {
	my ($seq,$choice) = @_;
	my ($Delta_H, $Delta_S) =(0,0);
	my $firstnucleotide = substr($seq,0,1);
	my $lastnucleotide = substr($seq,length($seq)-1,1);
		if($firstnucleotide eq 'G' || $firstnucleotide eq 'C') {
			$Delta_H += 0.1;
			$Delta_S += -2.8;
		} elsif($firstnucleotide eq 'A' || $firstnucleotide eq 'T') {
			$Delta_H += 2.3;
			$Delta_S += 4.1;
		} elsif(exists $IUPAC_DGN->{$firstnucleotide}){
			my $prob_gc = ($IUPAC_GC->{$firstnucleotide}/$IUPAC_DGN->{$firstnucleotide});
			my $prob_at = 1-$prob_gc;
				my $H_at = 2.3*$prob_at;
				my $S_at = 4.1*$prob_at;
				my $H_gc = 0.1*$prob_at;
				my $S_gc = -2.8*$prob_at;
			if($choice eq 'min'){
				my $min = $H_at < $H_gc ? $H_at : $H_gc;
				$Delta_H+=$min;
				$min = $S_at < $S_gc ? $S_at : $S_gc;
				$Delta_S+=$min;}
			else{
				my $max = $H_at > $H_gc ? $H_at : $H_gc;
				$Delta_H+=$max;
				$max = $S_at > $S_gc ? $S_at : $S_gc;
				$Delta_S+=$max;}
		}
		if($lastnucleotide eq 'G' || $lastnucleotide eq 'C') {
			$Delta_H += 0.1;
			$Delta_S += -2.8;
		} elsif($lastnucleotide eq 'A' || $lastnucleotide eq 'T') {
			$Delta_H += 2.3;
			$Delta_S += 4.1;
		} elsif(exists $IUPAC_DGN->{$lastnucleotide}){
			my $prob_gc = ($IUPAC_GC->{$lastnucleotide}/$IUPAC_DGN->{$lastnucleotide});
			my $prob_at = 1-$prob_gc;
				my $H_at = 2.3*$prob_at;
				my $S_at = 4.1*$prob_at;
				my $H_gc = 0.1*$prob_at;
				my $S_gc = -2.8*$prob_at;
			if($choice eq 'min'){
				my $min = $H_at < $H_gc ? $H_at : $H_gc;
				$Delta_H+=$min;
				$min = $S_at < $S_gc ? $S_at : $S_gc;
				$Delta_S+=$min;}
			else{
				my $max = $H_at > $H_gc ? $H_at : $H_gc;
				$Delta_H+=$max;
				$max = $S_at > $S_gc ? $S_at : $S_gc;
				$Delta_S+=$max;}
		}
	return (sprintf("%.2f",$Delta_H), sprintf("%.2f",$Delta_S));
}
#check for matches between passed comparing strings
#part of primer-dimer and hairpin sub routine
sub compare {
	my ($first,$second)= @_;
   	my ($match,$mismatch) = (0,0,0);
	my (@IUpacBase0,@IUpacBase1,$binding,$tmp);
	$binding='';
	my $pl = (length $first < length $second ? length $first : length $second);

   	for (0 .. $pl-1 ) {

		my $base0 = substr( $first, $_, 1);
		my $base1 = substr( $second, $_, 1);

		if(&is_base($base0) && &is_base($base1)) {
			if ($base0 ne $IUPAC_COMPL->{$base1}){
			       $tmp="0";
			} else{$tmp="1";$match++;}
		}elsif(&is_iupac_base($base0) && &is_base($base1)) {
			if($base0 eq 'M') {
				@IUpacBase0 = ('A','C');
			} elsif($base0 eq 'R') {
				@IUpacBase0 = ('A','G');
			} elsif($base0 eq 'W') {
				@IUpacBase0 = ('A','T');
			} elsif($base0 eq 'S') {
				@IUpacBase0 = ('C','G');
			} elsif($base0 eq 'Y') {
				@IUpacBase0 = ('C','T');
			} elsif($base0 eq 'K') {
				@IUpacBase0 = ('G','T');
			} elsif($base0 eq 'V') {
				@IUpacBase0 = ('A','C','G');
			} elsif($base0 eq 'H') {
				@IUpacBase0 = ('A','C','T');
			} elsif($base0 eq 'D') {
				@IUpacBase0 = ('A','G','T');
			} elsif($base0 eq 'B') {
				@IUpacBase0 = ('C','G','T');
			} elsif($base0 eq 'N') {
				@IUpacBase0 = ('A','C','G','T');
			}

			foreach my $base0 (@IUpacBase0) {
				if ($base0 ne $IUPAC_COMPL->{$base1}){
					$tmp="0";
				}else{$tmp="1";$match++;last;}
			}
		}elsif(&is_base($base0) && &is_iupac_base($base1)) {
			if($base1 eq 'M') {
				@IUpacBase1 = ('A','C');
			} elsif($base1 eq 'R') {
				@IUpacBase1 = ('A','G');
			} elsif($base1 eq 'W') {
				@IUpacBase1 = ('A','T');
			} elsif($base1 eq 'S') {
				@IUpacBase1 = ('C','G');
			} elsif($base1 eq 'Y') {
				@IUpacBase1 = ('C','T');
			} elsif($base1 eq 'K') {
				@IUpacBase1 = ('G','T');
			} elsif($base1 eq 'V') {
				@IUpacBase1 = ('A','C','G');
			} elsif($base1 eq 'H') {
				@IUpacBase1 = ('A','C','T');
			} elsif($base1 eq 'D') {
				@IUpacBase1 = ('A','G','T');
			} elsif($base1 eq 'B') {
				@IUpacBase1 = ('C','G','T');
			} elsif($base1 eq 'N') {
				@IUpacBase1 = ('A','C','G','T');
			}

			foreach my $base1 (@IUpacBase1) {
				if ($base1 ne $IUPAC_COMPL->{$base0}){
					$tmp="0";
				}else{$tmp="1";$match++;last;}
			}
		} elsif(&is_iupac_base($base0) && &is_iupac_base($base1)) {
			if($base0 eq 'M') {
				@IUpacBase0 = ('A','C');
			} elsif($base0 eq 'R') {
				@IUpacBase0 = ('A','G');
			} elsif($base0 eq 'W') {
				@IUpacBase0 = ('A','T');
			} elsif($base0 eq 'S') {
				@IUpacBase0 = ('C','G');
			} elsif($base0 eq 'Y') {
				@IUpacBase0 = ('C','T');
			} elsif($base0 eq 'K') {
				@IUpacBase0 = ('G','T');
			} elsif($base0 eq 'V') {
				@IUpacBase0 = ('A','C','G');
			} elsif($base0 eq 'H') {
				@IUpacBase0 = ('A','C','T');
			} elsif($base0 eq 'D') {
				@IUpacBase0 = ('A','G','T');
			} elsif($base0 eq 'B') {
				@IUpacBase0 = ('C','G','T');
			} elsif($base0 eq 'N') {
				@IUpacBase0 = ('A','C','G','T');
			}
			
			if($base1 eq 'M') {
				@IUpacBase1 = ('A','C');
			} elsif($base1 eq 'R') {
				@IUpacBase1 = ('A','G');
			} elsif($base1 eq 'W') {
				@IUpacBase1 = ('A','T');
			} elsif($base1 eq 'S') {
				@IUpacBase1 = ('C','G');
			} elsif($base1 eq 'Y') {
				@IUpacBase1 = ('C','T');
			} elsif($base1 eq 'K') {
				@IUpacBase1 = ('G','T');
			} elsif($base1 eq 'V') {
				@IUpacBase1 = ('A','C','G');
			} elsif($base1 eq 'H') {
				@IUpacBase1 = ('A','C','T');
			} elsif($base1 eq 'D') {
				@IUpacBase1 = ('A','G','T');
			} elsif($base1 eq 'B') {
				@IUpacBase1 = ('C','G','T');
			} elsif($base1 eq 'N') {
				@IUpacBase1 = ('A','C','G','T');
			}
			LOOP:foreach my $base0 (@IUpacBase0) {
				foreach my $base1 (@IUpacBase1) {
					if ($base0 ne $IUPAC_COMPL->{$base1}){
						$tmp="0";
					}else{$tmp="1";$match++;last LOOP;}
				}
			}
		}else{
			$tmp.="0";
		}
		$binding.=$tmp;
	}
	return ($match,$binding);
}
#finds the longest consecutive matches part of primer-dimer and hairpin sub-routine
#Arguments: binding primer sequence, binding pattern
sub find_longest_string {
	my ($primer,$binding)=@_;
	my ($lg_match,$string)='';
	my $start;
	my $new_binding = $binding;
	if ($binding !~ m/1111/){return (0,'','');}#should have less than 5 consecutive matched bases
	$lg_match = max map {length} $binding =~ /(1+)/g;
	my $tmp = '2' x $lg_match;
	$new_binding =~ s/1{$lg_match}/'2' x $lg_match/e;
	$start = index($new_binding,$tmp);
	$string = substr($primer,$start,$lg_match);
	return ($lg_match,$string,$new_binding);
}
#Computation of DeltaG kcal/mol using equation: dG = dH -TdS for the primer-dimer and hairpin formation
#Arguments: matched bases 
sub deltaG_primer_dimer {
	my $matched_bases = shift;
	my ($base0,$base,$H_min,$H_max,$S_min,$S_max,$salt_S,$dG_min,$dG_max,@temp,@IUpairVals_min,@IUpairVals_max);
	for(my $i=0; $i<length $matched_bases; $i++) {
		$base0 = substr($matched_bases,$i,1);
		$base = substr($matched_bases,$i+1,1);
		@temp = &calcIUpair($base0, $base, $i, 'min', $matched_bases);
		for(my $j=0; $j<3; $j++) {
			$IUpairVals_min[$j]+=$temp[$j];
		}
		@temp = &calcIUpair($base0, $base, $i, 'max', $matched_bases);
		for(my $j=0; $j<3; $j++) {
			$IUpairVals_max[$j]+=$temp[$j];
		}
	}
	$H_min = $IUpairVals_min[1];
	$S_min = $IUpairVals_min[2];
	my ($H_ends,$S_ends) = &HS_correct_ends($matched_bases,'min');
	$H_min+=$H_ends;
	$S_min+=$S_ends;
	$H_max = $IUpairVals_max[1];
	$S_max = $IUpairVals_max[2];
	($H_ends,$S_ends) = &HS_correct_ends($matched_bases,'max');
	$H_max+=$H_ends;
	$S_max+=$S_ends;
	$salt_S = &deltaS_correction($matched_bases);
	$S_min += $salt_S;
	$S_max += $salt_S;
	$dG_min = &recal_deltaG($H_min,$S_min);
	$dG_max = &recal_deltaG($H_max,$S_max);
	
	my $dG = ($dG_min < $dG_max ? $dG_min : $dG_max);

	return $dG;
}
#***************************************************************************************************************************
# Misc functions
#***************************************************************************************************************************
# Clean Blast directory- the individual files
# Arguments: none
sub clean_bl_files {
	system ("rm -r $TMP");
	return 0;
}
sub getIdFileData {
    my $idfile = shift;
    my %data;
    open(FILE,"<",$idfile) or die "ERROR: Could not open file $idfile: $! \n";
	while(<FILE>) {
	    chomp;
	    my @tmp = split(/\t/);
	    $data{$tmp[0]} = $tmp[1];
	}
    close(FILE);
    $data{filename} = &convertIntToString($data{filename}) if(exists $data{filename});
    return \%data;
}
sub get_siginfo_cds {
	my $ids = shift;
	my (@cols,@aalen,@nuclen,@start,@stop,@source,@all_function,@family,%tmp_function,$mean,@aalenMean);
	$mean = 0;
	open(CDS, "<", $CDS_TABLE) or die "ERROR: could not open cds table: $! \n";
	while(<CDS>) {
		chomp();
		@cols = split(/\t/);
		#($ref,$gi,$nclen,$aalen,$start,$end,$strand,$locus,$proid,$gorg,$phageclass)#order
		next unless(exists $ids->{$cols[0]});
		$ids->{$cols[0]} = [@cols[8,2,9,10,4,5,0]];
	}
	close(CDS);
	return;
}
sub getSigGenes {
	my ($file,$sig) = @_;
	my (@args,%genes,$good);
	$good = 0;
	open(FILE, "<", $file) or die "ERROR: could not open file $file: $! \n";
	while(<FILE>) {
		chomp();
		@args = split(/\t/);
		next unless(@args);
		if($args[0] eq $sig) {
			$good = 1;
		} elsif($args[0] =~ m/\w/ && $good) {
			last;
		}
		if($good) {
				$genes{$args[9]} = 1;
			}
	}
	close(FILE);
	return \%genes;
}
sub createSeqsFile {
	my ($genes,$file) = @_;
	my (@args);
	open(OUT, ">", $file) or die "ERROR: could not write to file: $! \n";
	open(IN, "<", $CDS_SEQS) or die "ERROR: could not open file: $! \n";
	while(<IN>) {
		chomp();
		@args = split(/\:/);
		if(exists $genes->{$args[0]}) {
			print OUT join("\n",'>'.$args[0],$args[1])."\n";
		}
	}
	close(IN);
	close(OUT);
}
sub condition_check{
	my $window = shift;
	my $condition;
	my @residues = split(//,$window);
	for (my	$i=0; $i < scalar @residues; $i++){
		my $base = $residues[$i];
		if(is_base($base)){
			$condition = 'deoxy';
		}
		if(is_iupac_base($base)){
			$condition = 'iupac';
			last;
		}
	}
	return $condition;
}
sub is_base {
	my $base = shift;
     	if (($base eq 'A') ||
		($base eq 'G') ||
		($base eq 'C') ||
		($base eq 'T')) {
		return 1;
	}
	return 0;
}
sub is_iupac_base {
	my $base = shift;
   	if (($base eq 'M') ||
		($base eq 'R') ||
		($base eq 'W') ||
		($base eq 'S') ||
		($base eq 'Y') ||
		($base eq 'K') ||
		($base eq 'V') ||
		($base eq 'H') ||
		($base eq 'D') ||
		($base eq 'B') ||
		($base eq 'N'))
	{
		return 1;
	}
	return 0;
}
sub count_AT {
	my $sequence = shift;
	my $count_AT = ($sequence =~ tr/A//);
	$count_AT += ($sequence =~ tr/T//);
	return $count_AT;
}
sub count_GC {
	my $sequence = shift;
	my $count_GC = ($sequence =~ tr/G//);
	$count_GC += ($sequence =~ tr/C//);
	return $count_GC;
}
sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
sub primer_min {
        my $primer = shift;
        $primer =~ s/[ATYRWKMDVHBN]/A/g;
        $primer =~ s/[CGS]/G/g;
	return $primer;
}
sub primer_max {
        my $primer = shift;
        $primer =~ s/[CGYRSKMDVHBN]/G/g;
        $primer =~ s/[ATW]/A/g;
	return $primer;
}
sub rev_complement {
	my $revrg_primer = shift;
	my $revcomp;
	my @residues = split(//,$revrg_primer);
	foreach my $r(@residues){
		$revcomp.=$IUPAC_COMPL->{$r};
	}
	return reverse($revcomp);#return 5'-3'
}
sub get_3Pri_5bases{
	my ($sequence) = @_;
	my $last5bases;
	$last5bases = substr($sequence,-5,5);#5'-3' last five bases
	return $last5bases;
}

