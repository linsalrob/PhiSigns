PhiSiGns
~~~~~~~~


DESCRIPTION
-----------
A computational tool for the identification of phage signature genes and PCR primer design


PURPOSE
-------
PhiSiGns is intended for the phage biologists who are interested  

1) in finding shared conserved genes in phages for evolution-based studies and/or

2) in designing primers on shared genes that can be used to discover unknown phages in the environment 



STANDALONE Vs. WEB-BASED VERSION
--------------------------------
PhiSiGns is a standalone and web-based application. Both versions provide the same features and results, but differ in the following:

1) The standalone version does not come with a phage genome database or precalculated BLASTP outputs. The user inputs the phage genomes to be compared for signature gene identification and primer design in an input txt file. The phage genomes are input as NCBI refseq numbers (see Standalone version - running mode). The GenBank files for the user-selected set of phage genomes are then obtained from NCBI and the gene annotations for the protein coding regions are imported accordingly. In contrast, the PhiSiGns web-based version consists of a database of 636 phages and 33 archaeal viruses (derived from the phage database on the PhAnToMe website, Feb 2011 [http://www.phantome.org/Downloads]), and pre-calculated BLASTP outputs for all these genomes (computed using PhiSiGns at an E-value cut-off of 10). The gene annotations are imported from the SEED [http://www.theseed.org/wiki/Home_of_the_SEED] and the proteins that lack annotation are extracted from GenBank. Therefore, with the standalone version users have the ability to compare phage genomes that are not part of the web-based PhiSiGns phage database.

2) The standalone version does not provide the option to select/deselect genes within the selected signature gene group for alignment and primer design. However, users can upload their own alignment containing only the sequences of interest.


STANDALONE VERSION - PREREQUISITES
----------------------------------
Standalone version is a perl script (phisigns.pl) and requires the following softwares and modules to run:

Resources
	Direct internet connection

Software
	BLASTALL (Blast standalone package)
	CLUSTALW 1.8 (Multiple sequence alignment tool)
	
BioPerl modules 
	Bio::Seq
	Bio::SeqIO
	Bio::SeqFeature::Generic
	Bio::DB::GenBank
	Bio::DB::Fasta
	Bio::Index::Fasta
	Bio::Tools::Run::StandAloneBlast
	Bio::Tools::Run::Alignment::Clustalw
	Bio::Align::Utilities qw(aa_to_dna_aln)
	Bio::AlignIO;

Perl modules
	Getopt::Long
	File::Basename
	List::Util qw(max)

Files
	input.txt

Also, you must have administrative rights to create, delete, read, or write files in the PhiSiGns working directory whereever located on your computer.
 

STANDALONE VERSION - SETTING VARIABLES
--------------------------------------
Following variables must be set before running phisigns.pl:

$BLASTALL (BLASTALL executable directory path)
Example: $BLASTALL  = '/home/user/blast/bin/blastall'

$FORMATDB (FORMATDB executable directory path)
Example: $FORMATDB = '/home/user/blast/bin/formatdb'

$ENV{CLUSTALDIR} (CLUSTALW directory path)
Example: $ENV{CLUSTALDIR} = '/home/user/clustalw1.8/' 


STANDALONE VERSION - COMMAND OPTIONS
------------------------------------

Option	  Description			   Default Argument
___________________________________________________________
-d	  PhiSiGns working directory path     -	   string
-m	  run mode			      0	      int
-e	  BLAST E-value cut-off		     10	    float
-c	  BLAST coverage cut-off	     10	    float
-i	  input file			      -	   string
-g	  selected SiG_# for primer design    -	   string
-a	  user uploaded alignment file	      -	   string
-aln_only Display default CLUSTALW alignment  -	     none
-id	  id of current job		      -	   string


STANDALONE VERSION - RUNNING MODE
---------------------------------
There are two modes to run PhiSiGns (phisigns.pl):


MODE 0
------

Identifies the signature genes amongst the user-selected phage genomes as listed in the input file ('-i') saved under the PhiSiGns working directory. The input file must include the NCBI refseq number for the phage genomes to be compared for signature gene identification and primer design. Type the following command line to execute phisigns.pl in Mode = 0:

Example:  % phisigns.pl -d /home/user/phisigns/ -m 0 -e 0.1 -c 25 -i input.txt -id test

BLAST E-value cut-off ('-e') and coverage ('-c') cut-off can be set by the user, if not specified a default value of 10 and 10 is used, respectively. Once the above command is executed and finished running, data and result files are generated with respect to the job id specified ('-id') in the command line under the PhiSiGns working directory.

Output files:
	data/.blastoutput			BLASTALL results for the selected phage genomes
	data/.cdstbl.ffd			Phage CDS table
	data/.cds_nucleotide.fa			CDS Nucleotide sequences
	output/.sig_gps.txt			short list of identified signature genes (SiG's)
	output/.sig_gps_detail.txt		detailed list of SiG's
	output/.log.txt				phisigns log file


MODE 1
------

Design PCR primer pairs for a selected signature gene (SiG_#) from the list of identified SiG's ('.sig_gps.txt' or .sig_gps_detail.txt' file) generated in Mode = 0. The primers are designed taking into account the minimum and maximum values specified for the primer parameters as in the pparams.txt file. Users can modify the values of these parameters to their preference. 

Primer Parameters options

Flag	 Description				Default		    Range
_________________________________________________________________________
minl	 min primer length			  16		    10-28
maxl	 max primer length			  28		    10-28
mingc	 min GC content				  30		    20-80
maxgc	 max GC content				  80		    30-80
minbtm	 min Basic Melting Temperature		  30		    30-80
maxbtm	 max Basic Melting Temperature		  80		    30-80
minstm	 min Salt Adjusted Temperature		  30		    30-80
maxstm	 max Salt Adjusted Temperature		  80		    30-80
minnnh	 min Nearest Neighbor Temperature	  30		    30-80
maxnnh	 max Nearest Neighbor Temperature	  80		    30-80
mindg	 min delta G for a primer		 -20		   -30-0
minpl	 min product length			 400		 100-2500
maxpl	 max product length			2000		 100-2500
dgn	 primer degeneracy 			1000		   1-1000
gcclamp	 3' GC clamp			           y		   y or n
max3stb  maximum 3' stability			   y		   y or n
compl	 Complementarity			   y		   y or n


Type the following command line to execute phisigns.pl in Mode = 1 to only display program generated alignment for a selected signature gene:

Example:  % phisigns.pl -d /home/user/phisigns/ -m 1 -g SiG_5 -aln_only -id test


Type the following command line to execute phisigns.pl in Mode = 1 with the program generated alignment:

Example:  % phisigns.pl -d /home/user/phisigns/ -m 1 -g SiG_5 -id test


Users can also upload their own nucleotide sequence alignment (in CLUSTALW format only with '.aln' file extension) for the selected signature gene group to design PCR primer pairs.  The user uploaded alignment file should be under the phisigns working directory. 

Type the following command line to execute phisigns.pl in Mode = 1 with a user-generated alignment:

Example:  % phisigns.pl -d /home/user/phisigns/ -m 1 -a test.aln -id test

Once the above command is executed and finished running, results files for the user uploaded alignment are generated with respect to the job id specified ('-id') in the command line under the phisigns working directory.

Output files:
	output/.seqs			sequences in FASTA format for the user-selected SiG_#
	output/.aln			nucleotide sequence CLUSTALW alignment for the selected SiG_#
	output/.pairs.txt		short list of potential primer pairs found for the selected SiG_#
	output/.pairs_detail.txt	detailed list of the potential primer pairs found for the selected SiG_#
	output/.log.txt			phisigns log file


Note: PhiSiGns (phisigns.pl) in Mode = 1 will execute only after PhiSiGns is executed in Mode = 0. Unless the selected set of phage genomes is different, Mode = 0 need not be executed more than once. To design PCR primer pairs with a user-uploaded alignment Mode = 1 can be executed without running Mode = 0 first. 


WEB VERSION -  RUNNING
----------------------
A. Identify signature genes
	1. Select and add phage genomes of interest from the list of available phage genomes
	2. Select BLAST E-value and coverage cut-off for similarity searches
	3. Click on "Get Signature Genes"
	4. View the list on the webpage or download the table on a local machine
B. Design PCR primer pairs
	1. Select a signature gene from the list of identified signature genes amongst the phage genomes of interest
	2. Select minimum and maximum primer parameter values [optional]
	3. Select/deselect genes within the selected signature gene group for alignment and primer design [optional]
	4. View program generated CLUSTALW alignment for selected genes [optional]
	5. Upload your own nucleotide sequence alignment for the selected signature gene [optional]
	6. Click on "Design Primers for Selected SiG Genes"
	7. View the results on the webpage or download them on a local machine

For more details on web version, please see http://phisigns.sourceforge.net/


WEB VERSION - DEPENDENCIES
--------------------------
The perl script require these Perl modules:
	CGI qw(:standard :html3)
	CGI::Carp qw(fatalsToBrowser)
	Fcntl qw(:DEFAULT :flock)


