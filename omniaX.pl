#!/usr/bin/perl -w
use Bio::SeqIO;
use Getopt::Long;

my $usage = "";

# define arguments and defaults
my ($forward, $reverse, $bam, $db);
my $circ = 1;
my $out = "";
my $spadesOpt = "";
GetOptions (
  "forward=s" => $forward, # if defined runs full pipeline
  "reverse=s" => $reverse, # optional, only can be defined with forward
  "spadesOpt=s" => \$spadesOpt, # string of spades options
  "circular!" => \$circ, # define if expected genomes can be circular
  "bam=s" => \$bam, # overwrites forward/reverse reads, analysis starts here
  "db=s" => \$db, # run analysis with known ME table
  "out=s" => \$out # output prefix
) or die("Error in command line arguments\n\n$usage\n");

my $inFile = shift;

my $in  = Bio::SeqIO->new(-file => $inFile, -format => 'Fasta');

while ( my $seqObj = $in->next_seq() ) {
  my $name = $seqObj->id;
  my $circlatorOut = "$name.circ.fa";
  my $getOrfOut = "$name.orfs.fa";
  open(OUT, "> $circlatorOut");
  print OUT ">", $name, "\n", $seqObj->seq, "\n";
  close(OUT);
  system("getorf -minsize 150 -circular Y -sequence $circlatorOut -outseq $getOrfOut");
}

#
# # load for metaspades
# module load spades
# # load for circlator
# module load Anaconda
# source activate myenv
# module load samtools
# module load mummer
# module load prodigal
# module load bwa
# # load for orf finder
# module load emboss
# # load for blastp
# module load blast+
# #load for hhsearch
# module load hhsearch
#
# # quick filter reads against human genome to reduce reads for de novo assembly
# # create initial contigs w/metaspades
# # complete circular genomes with circlator
# # find ORFs and translate
# getorf -minsize 150 -circular Y -sequence $circlatorOut -outseq $getOrfOut
# # blastp ORFs against curated viral db
#
# # input poorly matched or unmatched translated ORFs into hhsearch (hhsearch database should be restricted a curated set, PDB_mmCIF30)
#
# # last remaining unidentified ORFs should then be run into cdhit 90, 60, 30, 20? to annotate those that cluster (arbitrary id). If they don't cluster, toss 'em.
# # there should be some semblance of similarity even amongst unknowns, especially to have any chance of producing logical/accurate clusters
#
# # interpret blast, hhsearch, and cd-hit outputs, merge, and calculate mutual exclusivity/co-occurrence of protein categories
# ### output ORF hits by genome into table
# ### output ME table
# # build network clusters for ME table (statnet: betweenness, coreness, or cluster_fast_greedy)
# ### output cluster definitions
# # use these clusters to categorize input seqs (ie cluster contains which characteristics, calculate distance between seq and groups (set threshold for correct distance), include pre-defined groups for annotation)
# ### output hit table (txt and html, html will link to all other tables for easy navigation)
# ### visualization? Shiny?
# ##### if it all is stored in one table it becomes really easy to visualize

# output format(s)
#
# (1) merged ORF annotation table
## genome, start, end, orfid, annotation(s)
# (1.1) cd-hit clustering results
# (2) ME table
# (3) cluster definitions
## clusterid, list of proteins (score?)
# (4) final genome calls
# genome, clusterid, score

# proof of concept analysis. Load 10,100,...,n representatives of non-segmented genomes into this program. Calculate ROC AUC using TPR and FPR.
## TP = group is called correctly above confidence threshold
## FP = group is called incorrectly above confidence threshold
## FN = group is called correctly below confidence threshold
## TN = group is called incorrectly below confidence threshold
# most accurate learning set with be used for ME table annotation

# repeat analysis using best ME table with 1 excluded known

# repeat analysis on unknown dataset
