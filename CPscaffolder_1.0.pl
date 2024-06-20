#!/usr/bin/perl 
use strict;
###########################################################

=head1 SYNOPSIS

  Basic usage : perl CPscaffolder.pl -i contigs -r PacBioReads [-options]

  Example : 
       1. $ perl CPscaffolder.pl -i contigs.fasta -r pacbioreads.bas.h5
       2. $ perl CPscaffolder.pl -i contigs.fasta -b dir/blasr.m5 
       3. $ perl CPscaffolder.pl -i contigs.fasta -r pacbioreads.fastq -o scaffolds.fasta
       2. $ perl CPscaffolder.pl -i contigs.fasta -r pacbioreads.fastq -t 24 -g T 
       3. $ perl dir/CPscaffolder.pl -i dir/contigs.fa -r dir/pacbioreads.fq -c dir/config.ini 

  Option Description (default_value) :
       -h --help_message
            Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters.
       -v --version
            Print version number
       -i --input_contigs (contigs.fasta)
            Input fasta file containing contig sequences.
       -r --reference/pacbio_reads (reads.[fasta|fastq|bas.h5|pls.h5])
            Input pacbio reads file.
            bas.h5 and fastq formats are recommended.
       -o --output_file   (scaffolded.fasta)
            The final scaffolds, provided in FASTA format.
       -c --config_file   (software_directory/config.ini)
            Input configuration file, which is using the INI file format. 
       -t --threads (12)
            Number of threads (CPUs) to run BLASR and PGA.
       -b --blasr_m5_file
            Users could upload the alignment result in BLASR m5 format instead of PacBio reads.
       -g --gap_filling [T/F] (default : F)
            Using pacbio reads to fill gaps.

   Software Requirement :
       - BLASR (http://bix.ucsd.edu/projects/blasr/)
       - BioPerl (http://www.bioperl.org/wiki/Installing_BioPerl)
       - Statistics::Descriptive Perl module (http://www.cpan.org/modules/INSTALL.html)
       - Config::IniFiles Perl module (http://www.cpan.org/modules/INSTALL.html)

=head1 AUTHOR

  liqi , liqi@ihb.ac.cn

=head1 COPYRIGHT
  
  Copyright - Lab of Algal Genomics

=cut

###########################################################

use Bio::SeqIO;
use IO::String;
use threads;
use Config::IniFiles;
use Pod::Usage;
use Getopt::Long;
use File::Basename;
use Statistics::Descriptive;
use Cwd;

my ( $help, $contigs, $pacbio_reads, $outputfile, $version);
my $blasr_m5_file="0";
my $pacbio_reads="pacbio.fastq";
my $threads=24;
my $gap_on_off="F";
my $this_dir_path      = getcwd;
my $software_directory = dirname(__FILE__);
my $config             = "$software_directory/config.ini";

GetOptions(
    'h!'     => \$help,
    'v!'     => \$version,
    'i=s{1}' => \$contigs,
    'r=s{1}' => \$pacbio_reads,
    'c=s{1}' => \$config,
    'o=s{1}' => \$outputfile,
    'b=s{1}' => \$blasr_m5_file,
    't=i{1}' => \$threads,
    'g=s{1}' => \$gap_on_off,
);

pod2usage if $help;

if($version){
    print "\n\tCPscaffolder\nVersion : 1.0\nRelease date : 20-Aug-2014\n\n";
    exit;
}

(
    warn(
"\n[WARNING]\n\tMissing a required option!\n\tYou must appoint two input files (contigs and PacbioReads).\n\n\tRun with -h for a list of commands.\n\n"
      )
      and exit
) unless ( $contigs && $pacbio_reads || $contigs && $blasr_m5_file);

my $cfg = Config::IniFiles->new( -file => "$config" );

# Blasr Options:
my $minPctSimilarity    = $cfg->val( 'Blasr_Options', 'minPctSimilarity' );
my $minReadLength     = $cfg->val( 'Blasr_Options', 'minReadLength' );
my $minSubreadLength  = $cfg->val( 'Blasr_Options', 'minSubreadLength' );
my $maxScore          = $cfg->val( 'Blasr_Options', 'maxScore' );
my $bestn             = $cfg->val( 'Blasr_Options', 'bestn' );
my $nCandidates       = $cfg->val( 'Blasr_Options', 'nCandidates' );
my $minSeedLength     = $cfg->val( 'Blasr_Options', 'minSeedLength' );
my $advanceExactMatches=$cfg->val( 'Blasr_Options', 'advanceExactMatches' );
my $other_options     = $cfg->val( 'Blasr_Options', 'other_options' );
my $blasr_options = "--minPctSimilarity $minPctSimilarity --minReadLength $minReadLength --minSubreadLength $minSubreadLength --maxScore $maxScore --bestn $bestn --nCandidates $nCandidates --minMatch $minSeedLength --advanceExactMatches $advanceExactMatches $other_options";

# Alignments Screening Options:
my $min_contig_length = $cfg->val( 'Alignments_Screening', 'min_contig_length' );
my $min_mapQV         = $cfg->val( 'Alignments_Screening',   'min_mapQV' );
my $min_match_length  = $cfg->val( 'Alignments_Screening',   'min_match_length' );
my $min_FORSimilarity = $cfg->val( 'Alignments_Screening',   'min_FORSimilarity' );

# PGA Options:
my $pga_options =
    $cfg->val( 'PGA_Options', 'Iteration' ) . " "
  . $cfg->val( 'PGA_Options', 'RHO_value' ) . " "
  . $cfg->val( 'PGA_Options', 'BETA_value' ) . " "
  . $cfg->val( 'PGA_Options', 'Q0_value' );

# Scaffolding Options:
my $link_repetition_rate =
  $cfg->val( 'Scaffolding_Options', 'link_repetition_rate' );
my $candidate_paths_num =
  $cfg->val( 'Scaffolding_Options', 'candidate_paths_num' );
my $min_coverage = $cfg->val( 'Scaffolding_Options', 'min_coverage');

# Gap_filling Options:
my $min_gapfill_coverage = $cfg->val( 'GapFilling_Options', 'min_gapfill_coverage');
my $min_gapfill_FORSimilarity = $cfg->val( 'GapFilling_Options', 'min_gapfill_FORSimilarity');
my $max_gapfill_Reads = $cfg->val( 'GapFilling_Options', 'max_gapfill_Reads');

my @suffixlist = qw(.fasta .fastq .fa .fq .bas.h5 .pls.h5);
my ( $pacbio_name, $pacbio_path, $pacbio_suffix ) =
  fileparse( $pacbio_reads, @suffixlist );
my ( $contigs_name, $contigs_path, $contigs_suffix ) =
  fileparse( $contigs, @suffixlist );
my $pacbio_pathname = $pacbio_path . $pacbio_name;

my ($match_contigs_num,$unmatch_contigs_num,$totalNum_contigs,$totalNum_scaffolds,$work_dir,$tmpdir,$blasr_nproc,$pga_threads);
my (%contigs_id,%contigs_seq,%gap_fasta);
my (@matched_contigs,@matrix_score,@pairwise_contigs,@pairwise_ctgs_set,@distance_library,@orientation_library,@three_consecutive_link,@three_consecutive_link_lirary,@gap_info);

if($gap_on_off eq "T"){
    if($pacbio_suffix ne ".fastq" && $pacbio_suffix ne ".fq" || $blasr_m5_file ne "0"){
        print "\n[WARNING]\n\tCPscaffolder fills gap only by input a reference file (pacbio_reads) with fastq format.\n";
        print "\tNow resetting -g to F\n";
        $gap_on_off = "F";
    }
}


&function;
sub function{
    unless(-e "CPscaffolder_work_dir"){
        mkdir "CPscaffolder_work_dir" || die "can't create directory: PGAS_work_dir";
    }
    $work_dir = &gettime();
    my $date_time = &gettime();
    mkdir "CPscaffolder_work_dir/$work_dir" || die "can't create directory:PGAS_work_dir/$work_dir";
    mkdir "CPscaffolder_work_dir/$work_dir/TMPDIR" || die "can't create directory:TMPDIR";
    $tmpdir="CPscaffolder_work_dir/$work_dir/TMPDIR";

    unless($outputfile){
        $outputfile="$this_dir_path/CPscaffolder_work_dir/$work_dir/$contigs_name"."_scaffolds.fasta";
    }
    my ( $output_name, $output_path, $output_suffix ) = fileparse( $outputfile, @suffixlist );
    if($output_path eq "./"){
        $outputfile="$this_dir_path/CPscaffolder_work_dir/$work_dir/$output_name$contigs_suffix";
    }

    open (RECORD,">>CPscaffolder_work_dir/Jobs_info.txt");
    print RECORD "\n\n\n****************************************************\n";
    print RECORD "\t\t$date_time\n";
    print RECORD "****************************************************\n";
    print RECORD "$date_time\t[START]\n";

    open (LOG,">CPscaffolder_work_dir/$work_dir/log.txt");
    print LOG "\n\n\n****************************************************\n";
    print LOG "\t\t$date_time\n";
    print LOG "****************************************************\n";
    print LOG "$date_time\t[START]\n";

    print "\n[START] $date_time ****************************************************\n";
    print RECORD "$date_time\t[WorkDir]\t$work_dir\n";
    print LOG "$date_time\t[INPUT]\n\t\t\t\tContigs : $contigs_name$contigs_suffix\n\t\t\t\tPacbio Reads : $pacbio_name$pacbio_suffix\n\t\t\t\tConfig : $config\n";
    print RECORD "$date_time\t[INPUT]\n\t\t\t\tContigs : $contigs_name$contigs_suffix\n\t\t\t\tPacbio Reads : $pacbio_name$pacbio_suffix\n\t\t\t\tConfig : $config\n";
    print LOG "$date_time\t[CONFIGURATION]\n\t\t\t\tblasr options : $blasr_options\n\t\t\t\tminmum contig length : $min_contig_length bp\n\t\t\t\tminimum mapQV : $min_mapQV\n\t\t\t\tminimum match length : $min_match_length bp\n\t\t\t\tminimum FORSimilarity : $min_FORSimilarity\n\t\t\t\tminimum coverage : $min_coverage\n\t\t\t\tpga options : $pga_options\n\t\t\t\tlink repetition rate : $link_repetition_rate\n\t\t\t\tcandidate paths num : $candidate_paths_num\n\t\t\t\tnumber of threads : $threads\n";
    print RECORD "$date_time\t[CONFIGURATION]\n\t\t\t\tblasr options : $blasr_options\n\t\t\t\tminmum contig length : $min_contig_length bp\n\t\t\t\tminimum mapQV : $min_mapQV\n\t\t\t\tminimum match length : $min_match_length bp\n\t\t\t\tminimum FORSimilarity : $min_FORSimilarity\n\t\t\t\tminimum coverage : $min_coverage\n\t\t\t\tpga options : $pga_options\n\t\t\t\tlink repetition rate : $link_repetition_rate\n\t\t\t\tcandidate paths num : $candidate_paths_num\n\t\t\t\tnumber of threads : $threads\n";
    print "[CONFIGURATION]\n\tblasr options : $blasr_options\n\tminimum contig length : $min_contig_length bp\n\tminimum mapQV : $min_mapQV\n\tminimum match length : $min_match_length bp\n\tminimum FORSimilarity : $min_FORSimilarity\n\tminimum coverage : $min_coverage\n\tpga options : $pga_options\n\tlink repetition rate : $link_repetition_rate\n\tcandidate paths num : $candidate_paths_num\n\tnumber of threads : $threads\n";    
    print RECORD "$date_time\t[INFO] More details : $work_dir/log.txt\n";
    print RECORD "****************************************************\n\n";
    my $date_time = &gettime();
    print "[INFO] $date_time Now retrieving long contigs (>= $min_contig_length bp) from the target contigs ...\n";
    print LOG "$date_time\t[INFO] Now retrieving long contigs (>= $min_contig_length bp)\n";
    &contigs_screening;

    if($blasr_m5_file eq 0){
        my $date_time = &gettime();
        print "[INFO] $date_time Now running BLASR ...\n";
        print LOG "$date_time\t[INFO] Now running BLASR\n";
        &format if ($pacbio_suffix eq ".fasta" || $pacbio_suffix eq ".fa");
        &running_blasr;
    }else{
        print "[INFO] $date_time Now importing the BLASR.m5 file ...\n";
        print LOG "$date_time\t[INFO] Now importing the BLASR.m5 file\n";
        system("cut -d ' ' -f 1-17 $blasr_m5_file > $tmpdir/blasr.1.out");
        &blasr_primary_screening(1);
        unlink "$tmpdir/blasr.1.out";
        system("mv $tmpdir/blasr.1.pscreen $tmpdir/blasr.pscreen");
        # Add Title
        system("sed -i '1i qName\ttName\tqStrand\tqStart\tqEnd\tqLength\ttStrand\ttStart\ttEnd\ttLength\tnumMatch\tnumMismatch\tnumIns\tnumDel\tnumLeftMismatch\tnumRightMismatch\tmapQV\tscore\tFORSimilarity' $tmpdir/blasr.pscreen");
        # EndMark
        system("echo 'EndMark/0_1000\ttName\n' >> $tmpdir/blasr.pscreen");
    }

    my $date_time = &gettime();
    print "[INFO] $date_time Now screening alignment results ...\n";
    print LOG "$date_time\t[INFO] Now screening alignment results\n";
    &blasr_result_screening;

    &build_contigs_hash;

    my $date_time = &gettime();
    print "[INFO] $date_time Now calculating MeanDistance and Orientation ...\n";
    print LOG "$date_time\t[INFO] Now calculating MeanDistance and Orientation\n";
    &calculating_MeanDistance_and_Orientation;
    &build_three_consecutive_link_lirary;

    my $date_time = &gettime();
    print "[INFO] $date_time Now preparing to build Fitness Matrix ...\n";
    print LOG "$date_time\t[INFO] Now preparing to build Fitness Matrix\n";
    my $alignment = &check_alignments;
    my $date_time = &gettime();
    if(!$alignment){
        print "[ERROE] $date_time No alignments after screening!\n\n";
        print LOG "$date_time\t[ERROR] No alignments after screening!\n";
        exit;
    }

    my $date_time = &gettime();
    print "[INFO] $date_time Now building Fitness Matrix ...\n";
    print LOG "$date_time\t[INFO] Now building Fitness Matrix\n";
    &bulid_fitness_matrix;

    my $date_time = &gettime();
    print "[INFO] $date_time Now running PGA ...\n";
    print LOG "$date_time\t[INFO] Now running PGA\n";
    &running_PGA;

    my $date_time = &gettime();
    print "[INFO] $date_time Now evaluating the link path ...\n";
    print LOG "$date_time\t[INFO] Now evaluating the link path\n";
    &connection_paths_evaluate;

    if($gap_on_off eq "T"){
        my $date_time = &gettime();
        print "[INFO] $date_time Now generating gaps sequence ...\n";
        print LOG "$date_time\t[INFO] Now generating gaps sequence\n";
        &get_gap_fasta;
    }

    my $date_time = &gettime();
    print "[INFO] $date_time Now generating scaffolds sequence ...\n";
    print LOG "$date_time\t[INFO] Now generating scaffolds sequence\n";
    &congtigs_scaffolding;

    my $date_time = &gettime();
    print "[INFO] $date_time The output file is : $outputfile\n";
    print "[INFO] $date_time \n";
    my @pre_statistics=&statistics("$tmpdir/$contigs_name.newseq.fasta");
    my @post_statistics=&statistics("$outputfile");
    my $total_contigs=sprintf("%12d",$totalNum_contigs);
    my $total_scaffolds=sprintf("%12d",$totalNum_scaffolds);

    open(REPORT,">CPscaffolder_work_dir/$work_dir/report.table");
    print"\t\t\t\t   CPscaffolder Report\n";
    print REPORT "\t\t\t   CPscaffolder Report\n";
    print "\t\t---------------------------------------------------\n";
    print REPORT "\t---------------------------------------------------\n";
    print "\t\t\t\t|\tpre\t|\tpost\n";
    print REPORT "\t\t\t|\tpre\t|\tpost\n";
    print "\t\t---------------------------------------------------\n";
    print REPORT "\t---------------------------------------------------\n";
    if($unmatch_contigs_num > 0){
        print "\t\tNO. contigs\t| $total_contigs\t| $total_scaffolds *\n";
        print REPORT "\tNO. contigs\t| $total_contigs\t| $total_scaffolds *\n";
    }else{
        print "\t\tNO. contigs\t| $total_contigs\t| $total_scaffolds\n";
        print REPORT "\tNO. contigs\t| $total_contigs\t| $total_scaffolds\n";
    }
    print "\t\tTotal length\t| $pre_statistics[0]\t| $post_statistics[0]\n";
    print REPORT "\tTotal length\t| $pre_statistics[0]\t| $post_statistics[0]\n";
    print "\t\tLargest contig\t| $pre_statistics[1]\t| $post_statistics[1]\n";
    print REPORT "\tLargest contig\t| $pre_statistics[1]\t| $post_statistics[1]\n";
    print "\t\tAverage length\t| $pre_statistics[2]\t| $post_statistics[2]\n";
    print REPORT "\tAverage length\t| $pre_statistics[2]\t| $post_statistics[2]\n";
    print "\t\tGC (%)\t\t| $pre_statistics[3]\t| $post_statistics[3]\n";
    print REPORT "\tGC (%)\t\t| $pre_statistics[3]\t| $post_statistics[3]\n";
    print "\t\tN50\t\t| $pre_statistics[4]\t| $post_statistics[4]\n";
    print REPORT "\tN50\t\t| $pre_statistics[4]\t| $post_statistics[4]\n";
    # print "\t\tN90\t\t| $pre_statistics[5]\t| $post_statistics[5]\n";
    # print REPORT "\tN90\t\t| $pre_statistics[5]\t| $post_statistics[5]\n";
    print "\t\t---------------------------------------------------\n"; 
    print REPORT "\t---------------------------------------------------\n";
    print "\t\tAll statistics are based on contigs of size >= $min_contig_length bp.\n";
    print REPORT "\tAll statistics are based on contigs of size >= $min_contig_length bp.\n";
    print "\t\t * NOTE: Post NO. Contigs included $unmatch_contigs_num unmatched contigs.\n" if ($unmatch_contigs_num > 0);
    print REPORT "\t * NOTE: Post NO. Contigs included $unmatch_contigs_num unmatched contigs.\n" if($unmatch_contigs_num > 0);

    close REPORT;

    print LOG "$date_time\t[OUTPUT]\n\t\t\t\t$outputfile\n";

    print LOG "$date_time\t[RESULT]\n";
    print LOG "\t\t\t\tNO. Scaffolds : $totalNum_scaffolds (include $unmatch_contigs_num unmatched contigs)\n";
    print LOG "\t\t\t\tNO. Original Contigs : $totalNum_contigs (longer than $min_contig_length bp)\n";

    print LOG "$date_time\t[INFO] The program run successfully !\n";

    print "\n[END] $date_time ****************************************************\n\n";
    print LOG "$date_time\t[END]\n";
    print LOG "****************************************************\n\n";

    close LOG;
    close RECORD;
}

# Select longer sequences ( >= min_contig_length bp) from the input contigs;
sub contigs_screening {
    my $i   = 0;
    my $in  = Bio::SeqIO->new( -file => "$contigs", -format => 'fasta' );
    open(NewSeq,">$tmpdir/$contigs_name.newseq.fasta");

    while ( my $seq = $in->next_seq() ) {
        my $length = $seq->length;
        if ( $length >= $min_contig_length ) {
            my $id = $seq->id;
            $id=~s/\|/_/g;
            $id=~s/_$//;
            my $sequence = $seq->seq;
            print NewSeq ">$id\n$sequence\n";
            $i++;
        }
    }

    $totalNum_contigs = $i;
    close NewSeq;
}

sub format {
    my $in  = Bio::SeqIO->new( -file => "$pacbio_reads", -format => 'fasta' );
    open(Newseq,">$pacbio_name.newseq.fasta");
    while ( my $seq = $in->next_seq() ) {
        my $id = $seq->id;
        my $seqence = $seq->seq;
        print Newseq ">$id\n$seqence\n";
    }
    close Newseq;
    rename $pacbio_reads, "$pacbio_reads.bak";
    rename "$pacbio_name.newseq.fasta", $pacbio_reads;
}

sub gettime {
    my ( $sec, $min, $hour, $mday, $mon, $year ) = localtime(time);
    ( $sec, $min, $hour, $mday, $mon, $year ) = (
        sprintf( "%02d", $sec ),
        sprintf( "%02d", $min ),
        sprintf( "%02d", $hour ),
        sprintf( "%02d", $mday ),
        sprintf( "%02d", $mon + 1 ),
        $year + 1900 );
    my $now_date = join( "-", ( $year, $mon, $mday ) );
    my $now_time = join( ":", ( $hour, $min, $sec ) );
    my $date_time = join( "T", ( $now_date, $now_time ) );
    return $date_time;
}

sub running_blasr {
    $blasr_nproc=int($threads / 3);
    my $remainder=$threads % 3;
    $blasr_nproc++ if($remainder!=0);

    if ( $pacbio_suffix eq ".bas.h5" || $pacbio_suffix eq ".pls.h5" ) {
        &blasr;
    }else {
        my ($totalline)= `wc -l $pacbio_reads`=~/^(\d+)/;
        my $linenum = (int($totalline/12) + 1)*4;
        system("split -l $linenum -a 1 -d $pacbio_reads $pacbio_name.");
        system("mv $pacbio_name.0 $tmpdir/$pacbio_name.1$pacbio_suffix");
        system("mv $pacbio_name.1 $tmpdir/$pacbio_name.2$pacbio_suffix");
        system("mv $pacbio_name.2 $tmpdir/$pacbio_name.3$pacbio_suffix");
        &blasr;
        unlink glob "$tmpdir/$pacbio_name.?$pacbio_suffix";
    }
    system("cat $tmpdir/blasr.1.pscreen $tmpdir/blasr.2.pscreen $tmpdir/blasr.3.pscreen > $tmpdir/blasr.pscreen");
    unlink glob "$tmpdir/blasr.?.pscreen";

    # Add Title
    system("sed -i '1i qName\ttName\tqStrand\tqStart\tqEnd\tqLength\ttStrand\ttStart\ttEnd\ttLength\tnumMatch\tnumMismatch\tnumIns\tnumDel\tnumLeftMismatch\tnumRightMismatch\tmapQV\tscore\tFORSimilarity' $tmpdir/blasr.pscreen");
    # EndMark
    system("echo 'EndMark/0_1000\ttName\n' >> $tmpdir/blasr.pscreen");
}

# Multithreading optimization for sequence alignment
sub blasr {
    my $thr1 = threads->new( \&blasr_command, 1 );
    my $thr2 = threads->new( \&blasr_command, 2 );
    my $thr3 = threads->new( \&blasr_command, 3 );
    $thr1->join;
    $thr2->join;
    $thr3->join;
}

sub blasr_command {
    my ($i) = @_;
    if ( $pacbio_suffix eq ".bas.h5"){
        system("blasr $pacbio_pathname.$i.bax.h5 $tmpdir/$contigs_name.newseq.fasta $blasr_options --nproc $blasr_nproc -m 5 --out $tmpdir/blasr.$i.m5");
    }else{
        system("blasr $tmpdir/$pacbio_name.$i$pacbio_suffix $tmpdir/$contigs_name.newseq.fasta $blasr_options --nproc $blasr_nproc -m 5 --out $tmpdir/blasr.$i.m5");
    }
    system("cut -d ' ' -f 1-17 $tmpdir/blasr.$i.m5 > $tmpdir/blasr.$i.out");
    unlink "$tmpdir/blasr.$i.m5";
    &blasr_primary_screening($i);
    unlink "$tmpdir/blasr.$i.out";
}

sub blasr_primary_screening{
    my ($i) = @_;
    open( BlasrResult, "$tmpdir/blasr.$i.out" ) || die "[Error]\n\tCannot find the file: ./$tmpdir/blasr.$i.out\n";
    open( Pscreening, ">$tmpdir/blasr.$i.pscreen" ) || die;
    while (<BlasrResult>) {
        chomp;
        my (
            $qName,       $qLength, $qStart,  $qEnd,
            $qStrand,     $tName,   $tLength, $tStart,
            $tEnd,        $tStrand, $score,   $numMatch,
            $numMismatch, $numIns,  $numDel,  $mapQV
        ) = split /\s+/;

        die "No alignments !\n" if (!$qName && $.==1);

        # Skip alignments that mapQV is below given parameters by the user
        next if ( $mapQV < $min_mapQV );

        if ( $pacbio_suffix eq ".bas.h5" || $pacbio_suffix eq ".pls.h5" ) {
            my ( $rStart, $rEnd ) = $qName =~ /(\d+)_(\d+)$/;
            $qStart  = $qStart - $rStart;
            $qEnd    = $qEnd - $rStart;
            $qLength = $rEnd - $rStart;
        }

        # Skip alignments that match length is below given parameters by the user
        next if ( $numMatch < $min_match_length );

        if ( $tStrand eq "-" ) {
            my $temp = $tStart;
            $tStart = $tLength - $tEnd;
            $tEnd   = $tLength - $temp;
        }

        # Calculate FORSimilarity
        my $numLeftMismatch = $qStart < $tStart ? $qStart : $tStart;
        my $numRightMismatch = ( $qLength - $qEnd ) < ( $tLength - $tEnd ) ? ( $qLength - $qEnd ) : ( $tLength - $tEnd );
        my $FORSimilarity = ($numMatch / ($numMatch +$numMismatch + $numIns + $numDel +$numLeftMismatch +$numRightMismatch)) * 100;

        # Skip alignments that FORSimilarity is below given parameter by the user
        next if ( $FORSimilarity < $min_FORSimilarity );

        $FORSimilarity = sprintf( "%.2f", $FORSimilarity );
        my $line = "$qName\t$tName\t$qStrand\t$qStart\t$qEnd\t$qLength\t$tStrand\t$tStart\t$tEnd\t$tLength\t$numMatch\t$numMismatch\t$numIns\t$numDel\t$numLeftMismatch\t$numRightMismatch\t$mapQV\t$score\t$FORSimilarity";
        print Pscreening "$line\n";
    }

    close BlasrResult;
    close Pscreening;
}

sub blasr_result_screening {
    open( Blasrpscreen, "$tmpdir/blasr.pscreen" ) || die "[Error]\n\tCannot find the file: ./$tmpdir/blasr.pscreen\n";
    open( BlasrScreen, ">$tmpdir/blasr.screen" ) || die "[Error]\n\tCannot find the file: ./$tmpdir/blasr.screen\n";
    print BlasrScreen "qName\ttName\tqStrand\tqStart\tqEnd\tqLength\ttStrand\ttStart\ttEnd\ttLength\tnumMatch\tnumMismatch\tnumIns\tnumDel\tnumLeftMismatch\tnumRightMismatch\tmapQV\tscore\tFORSimilarity\n";
    my $last;
    my @array;

    while (<Blasrpscreen>) {
        next if $.==1;
        chomp;
        my ($qName,$tName)=$_=~/^(\S+)\s+(\S+)/;
        if ($. == 2){
            $last = $qName;
            die "No alignments in blasr.pscreen !\n" if (!$qName);
        }

        if ( $last ne $qName ) {
            &pairwise_alignment_screening(@array) if $#array > 0;
            @array = ();
        }
        push( @array, $_ ); 
        $last = $qName;
    }

    close Blasrpscreen;
    close BlasrScreen;
}

sub pairwise_alignment_screening {
    my (@array) = @_;
    my @screen;
    foreach (@array) {
        my $lable = 0;
        my (
            $qNameA,            $tNameA,    $qStrandA,
            $qStartA,           $qEndA,     $qLengthA,
            $tStrandA,          $tStartA,   $tEndA,
            $tLengthA,          $numMatchA, $numMismatchA,
            $numInsA,           $numDelA,   $numLeftMismatchA,
            $numRightMismatchA, $mapQVA,    $scoreA,
            $FORSimilarityA
        ) = split /\s+/;

        foreach (@array) {
            my (
                $qNameB,            $tNameB,    $qStrandB,
                $qStartB,           $qEndB,     $qLengthB,
                $tStrandB,          $tStartB,   $tEndB,
                $tLengthB,          $numMatchB, $numMismatchB,
                $numInsB,           $numDelB,   $numLeftMismatchB,
                $numRightMismatchB, $mapQVB,    $scoreB,
                $FORSimilarityB
            ) = split /\s+/;
            if ( $tNameA eq $tNameB ) {
                $lable = 1 if ( $FORSimilarityA < $FORSimilarityB );
            }else {
                if (   $qStartA >= $qStartB && $qStartA <= $qEndB || $qEndA <= $qEndB && $qEndA >= $qStartB ){
                    my @sorted = sort {$a <=> $b } ( $qStartA, $qStartB, $qEndA, $qEndB );
                    my $match_overlap_length = $sorted[2] - $sorted[1];
                    $lable = 1 if ( $match_overlap_length > 500 && $scoreA > $scoreB ); 
                    $lable = 1 if ( $match_overlap_length >= 0.9 * $tLengthA);
                }
            }
        }

        push( @screen, $_ ) if ( $lable == 0 );
    }
    &calculating_distance(@screen) if ( @screen > 1 );
    &build_three_consecutive_link(@screen) if ( @screen > 2 );
}

sub calculating_distance {
    my (@array) = @_;
    my @sorted_array = map {$_->[0] }sort { $a->[1] <=> $b->[1] } map {[ $_, ( split /\s+/ )[3] ] } @array;
    my @comparison = @sorted_array;

    foreach(@sorted_array){
        print BlasrScreen "$_\n";
        my (
            $qNameA,            $tNameA,    $qStrandA,
            $qStartA,           $qEndA,     $qLengthA,
            $tStrandA,          $tStartA,   $tEndA,
            $tLengthA,          $numMatchA, $numMismatchA,
            $numInsA,           $numDelA,   $numLeftMismatchA,
            $numRightMismatchA, $mapQVA,    $scoreA,
            $FORSimilarityA
        ) = split /\s+/;
        push @matched_contigs,$tNameA if((grep{/$tNameA$/}@matched_contigs)==0);

        shift @comparison;
        foreach (@comparison) {
            my (
            $qNameB,            $tNameB,    $qStrandB,
            $qStartB,           $qEndB,     $qLengthB,
            $tStrandB,          $tStartB,   $tEndB,
            $tLengthB,          $numMatchB, $numMismatchB,
            $numInsB,           $numDelB,   $numLeftMismatchB,
            $numRightMismatchB, $mapQVB,    $scoreB,
            $FORSimilarityB
            ) = split /\s+/;

            my $adjacency_distance = ( $qStartB - $numLeftMismatchB ) - ( $qEndA + $numRightMismatchA );
            next if ( $adjacency_distance < -300 );
            my $start = $qEndA + $numRightMismatchA;
            my $end   = $qStartB - $numLeftMismatchB;

            # calculate insertion and deletion rate
            my $total_numMatch = $numMatchA + $numMismatchA + $numInsA + $numDelA + $numMatchB + $numMismatchB + $numInsB + $numDelB;
            my $Ins_rite = ( $numInsA + $numInsB ) / $total_numMatch;
            my $Del_rite = ( $numDelA + $numDelB ) / $total_numMatch;
        
            # calculate corrective distance
            my $corrective_distance = int( $adjacency_distance * (1 - $Ins_rite + $Del_rite));

            # calculate Connectivity Quality (CQ) Value
            my $error_probabilities = 1 - ($FORSimilarityA / 100) * ($FORSimilarityB / 100);
            # Can't take log of 0
            if($error_probabilities==0){
                $error_probabilities = 1 - 0.99*0.99;
            }
            my $CQ = -10 * log($error_probabilities) / log(10);
            $CQ = sprintf( "%.2f", $CQ );

            push @gap_info,"$qNameB\t$start\t$end\t$tNameA\t$tNameB\t$CQ" if ($FORSimilarityA >= $min_gapfill_FORSimilarity && $FORSimilarityB >= $min_gapfill_FORSimilarity);

            my $orientation = $tStrandA . $tStrandB;
            my $orientation_code;
            if ( $orientation eq "++" ) {
                $orientation_code = 1;
            }elsif ( $orientation eq "+-" ) {
                $orientation_code = 2;
            }elsif ( $orientation eq "-+" ) {
                $orientation_code = 3;
            }elsif ( $orientation eq "--" ) {
                $orientation_code = 4;
            }

            push @pairwise_ctgs_set,"$tNameA\t$tNameB\t$corrective_distance\t$orientation_code\t$CQ";
            push @pairwise_contigs, "$tNameA\t$tNameB" if ( ( grep {/$tNameB/} grep {/$tNameA/} @pairwise_contigs ) == 0 );
        }
    }
    print BlasrScreen "\n";
}

sub build_three_consecutive_link{
    my (@array) = @_;
    my @sorted_array = map { $_->[0] }sort { $a->[1] <=> $b->[1] } map { [ $_, ( split /\s+/ )[3] ] } @array;
    my $lable =1;
    my $long_link;

    for ( my $i = 0 ; $i < $#sorted_array ; $i++ ) {
        my (
            $qNameA,            $tNameA,    $qStrandA,
            $qStartA,           $qEndA,     $qLengthA,
            $tStrandA,          $tStartA,   $tEndA,
            $tLengthA,          $numMatchA, $numMismatchA,
            $numInsA,           $numDelA,   $numLeftMismatchA,
            $numRightMismatchA, $mapQVA,    $scoreA,
            $FORSimilarityA
        ) = split /\s+/, $sorted_array[$i];
        my (
            $qNameB,            $tNameB,    $qStrandB,
            $qStartB,           $qEndB,     $qLengthB,
            $tStrandB,          $tStartB,   $tEndB,
            $tLengthB,          $numMatchB, $numMismatchB,
            $numInsB,           $numDelB,   $numLeftMismatchB,
            $numRightMismatchB, $mapQVB,    $scoreB,
            $FORSimilarityB
        ) = split /\s+/, $sorted_array[ $i + 1 ];

        my $adjacency_distance = ( $qStartB - $numLeftMismatchB ) - ( $qEndA + $numRightMismatchA );
        next if ( $adjacency_distance < -300 );
        
        if($lable == 1){
            $long_link = $tNameA ."->". $tNameB;
        }else{
            $long_link = $long_link ."->". $tNameB;
        }
        $lable ++;
    }

    if ($lable==3){
        push @three_consecutive_link,$long_link;
    }elsif($lable>3){
        my @link = split /->/,$long_link;
        for(my $j=0;$j<$#link-1;$j++){
            my $threepair = $link[$j] ."->".$link[$j+1] ."->". $link[$j+2];
            push @three_consecutive_link,$threepair;
        }
    }
}

sub build_contigs_hash{
    my $in = Bio::SeqIO->new( -file => "$tmpdir/$contigs_name.newseq.fasta", -format => 'fasta' );
    my $i  = 1;
    my $j  = 0;
    while ( my $seq = $in->next_seq() ) {
        my $lable=0;
        foreach(@matched_contigs){
            if($_ eq $seq->id){
                $contigs_seq{ $seq->id } = $seq->seq;
                $i = sprintf( "%05d", $i );
                $contigs_id{$seq->id}= $i;
                $i++;
                $lable=1;
                last;
            }
        }
        if($lable==0){
            my $out = Bio::SeqIO->new(-file   => ">>$tmpdir/PacBio.unmatched.fasta",-format => 'fasta');
            $out->write_seq($seq);
            $j++;
        }
    }
    $match_contigs_num   = $i -1;
    $unmatch_contigs_num = $j;
}

sub calculating_MeanDistance_and_Orientation {
    foreach (@pairwise_contigs) {
        my ( $contigA, $contigB ) = $_ =~ /(\S+)\s+(\S+)/;
        
        my @pairwise_contigAB = grep { /$contigB/ } grep { /$contigA/ } @pairwise_ctgs_set;
        my @distanceAB = map { $_->[1] } map { [ $_, ( split /\s+/ )[2] ] } @pairwise_contigAB;

        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data( \@distanceAB );
        my $Q1                 = $stat->quantile(1);
        my $Q3                 = $stat->quantile(3);
        my $lower_inner_fences = $Q1 - 1.5 * ( $Q3 - $Q1 );
        my $upper_inner_fences = $Q3 + 1.5 * ( $Q3 - $Q1 );

        my ( @filtered_distance, @orientation_AB, @orientation_BA );
        foreach (@pairwise_contigAB) {
            my ( $contig_i, $contig_j, $distance, $orientation_code, $CQ ) = split /\s+/;
            if ($distance >= $lower_inner_fences && $distance <= $upper_inner_fences ){
                push @filtered_distance, $distance;

                my $i = $contigs_id{$contig_i} -1;
                my $j = $contigs_id{$contig_j} -1;
                $matrix_score[$i][$j] = $CQ + $matrix_score[$i][$j] + 2;
                $matrix_score[$j][$i] = $CQ + $matrix_score[$i][$j];

                if ( $contigA eq $contig_i ) {
                    push @orientation_AB, $orientation_code;
                }else {
                    push @orientation_BA, $orientation_code;
                }
            }
        }

        my $stat1 = Statistics::Descriptive::Full->new();
        $stat1->add_data( \@filtered_distance );
        my $mean_distance = int( $stat1->mean() );
        push @distance_library, "$contigA\t$contigB\t$mean_distance";
        
        if ( @orientation_AB > 0 ) {
            my $orientationAB_code = &orient_by_majority_rule(@orientation_AB);
            push @orientation_library,"$contigA\t$contigB\t$orientationAB_code";
        }
        if ( @orientation_BA > 0 ) {
            my $orientationBA_code = &orient_by_majority_rule(@orientation_BA);
            push @orientation_library,"$contigB\t$contigA\t$orientationBA_code";
        }
        &check_orientation_library( $contigA,$contigB );
    }
}

sub check_orientation_library{
    my ($contigi,$contigj) = @_;
    my @array = grep { /$contigi/ } grep { /$contigj/ } @orientation_library;
    if(@array==1){
        my($contigA,$contigB,$orientationAB_code)=split /\t/,$array[0];
        my $orientationBA_code;
        if($orientationAB_code eq "F_F"){
            $orientationBA_code = "R_R";
        }elsif($orientationAB_code eq "F_R"){
            $orientationBA_code = "F_R";
        }elsif($orientationAB_code eq "R_F"){
            $orientationBA_code = "R_F";
        }elsif($orientationAB_code eq "R_R"){
            $orientationBA_code = "F_F";
        }
        push @orientation_library, "$contigB\t$contigA\t$orientationBA_code";
    }
}

sub orient_by_majority_rule{
    my (@array) = @_;
    my ( $mode, $orientation_code );
    if ( @array == 1 ) {
        $mode = $array[0];
    }else {
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data( \@array );
        $mode = $stat->mode();
    }
    if ( $mode == 1 ) {
        $orientation_code = "F_F";
    }elsif ( $mode == 2 ) {
        $orientation_code = "F_R";
    }elsif ( $mode == 3 ) {
        $orientation_code = "R_F";
    }elsif ( $mode == 4 ) {
        $orientation_code = "R_R";
    }
    return $orientation_code;
}

sub check_alignments {
    open( BlasrScreen, "$tmpdir/blasr.screen" ) || die "[Error]\n\tCannot find the file: ./$tmpdir/blasr.screen\n";
    readline BlasrScreen;
    my $line=<BlasrScreen>;
    return ($line);
    close BlasrScreen;
}

sub bulid_fitness_matrix {
    &fitness_matrix_supplement;
    my $max = 0;
    for ( my $i = 0 ; $i < $match_contigs_num ; $i++ ) {
        for ( my $j = 0 ; $j < $match_contigs_num ; $j++ ) {
            $max = $matrix_score[$i][$j] if ( $matrix_score[$i][$j] > $max );
        }
    }
    $max = $max + 1;

    open( FitnessMatrix, ">$tmpdir/Fitness.matrix" );
    print FitnessMatrix "$match_contigs_num\n";
    
    my $min_CQ = -10 * log(1-($min_FORSimilarity/100)*($min_FORSimilarity/100)) / log(10);
    my $min_matrix_score = $min_coverage * $min_CQ;
    $min_matrix_score = sprintf( "%.2f", $min_matrix_score);

    for ( my $i = 0 ; $i < $match_contigs_num ; $i++ ) {
        for ( my $j = 0 ; $j < $match_contigs_num ; $j++ ) {
            $matrix_score[$i][$j] = 0 if($matrix_score[$i][$j] <= $min_matrix_score);
            $matrix_score[$i][$j] = $max - $matrix_score[$i][$j];
            $matrix_score[$i][$j] = sprintf( "%.0f", $matrix_score[$i][$j] );
            print FitnessMatrix "$matrix_score[$i][$j]\t";
        }
        print FitnessMatrix "\n";
    }
    close FitnessMatrix;
    @matrix_score=();
}

sub fitness_matrix_supplement{
    my @repetitive_contigs;
    foreach(@three_consecutive_link_lirary){
        my $line = $_;
        my ($left,$center,$right)=split /->/,$line;
        if((grep{/$left/}grep{/$center/}grep{/$right/}@three_consecutive_link_lirary)==1){
            push @three_consecutive_link_lirary,$right."->$center->".$left;
        }
        if((grep{/->$center->/}@three_consecutive_link_lirary) > 2){
            if((grep{/$center/}@repetitive_contigs) == 0){
                my $i = $contigs_id{$left} - 1;
                my $j = $contigs_id{$center} - 1;
                my $k = $contigs_id{$right} - 1;
                # 50 : -10log(1-0.99*0.99)/log(10) * 3 coverage
                $matrix_score[$i][$j] += 50;
                $matrix_score[$j][$i] += 50;
                $matrix_score[$j][$k] += 50;
                $matrix_score[$k][$j] += 50;
                push @repetitive_contigs,$line;
            }else{
                if((grep{/$left/}@repetitive_contigs) == 0){
                    my $i = $contigs_id{$left} - 1;
                    my $k = $contigs_id{$right} - 1;
                    $matrix_score[$i][$k] += 50;
                    $matrix_score[$k][$i] += 50;
                }
            }
        }
    }
}

sub running_PGA {
    my $max_ctgs_num = $match_contigs_num + 1;
    system("sed 's/const int cities = 600/const int cities = $max_ctgs_num/' $software_directory/pga.cpp > $tmpdir/pga_1.cpp");
    system("g++ -o $tmpdir/pga $tmpdir/pga_1.cpp");
    unlink "$tmpdir/pga_1.cpp";

    my @thread_array;
    my $pga_out;
    $pga_threads = $threads;
    if ($pga_threads > $candidate_paths_num/2){
        $pga_threads = int($candidate_paths_num/2);
    }
    for(my $t=0;$t<$pga_threads;$t++){
        my $x=int($candidate_paths_num / $pga_threads);
        $x=$x + ($candidate_paths_num % $pga_threads) if $t==0;
        $thread_array[$t] = threads->new( \&PGA_command,$t,$x );
        sleep 1;
        $pga_out .= "$tmpdir/pga.out.$t ";
    }
    foreach my $thread(@thread_array){
        $thread->join();
    }
    system("cat $pga_out > $tmpdir/ConnectionPathSets.txt");
    unlink glob "$tmpdir/pga.out.*";
}

sub PGA_command {
    my ($i,$j) = @_;
    system("$tmpdir/pga $tmpdir/Fitness.matrix $tmpdir/pga.out.$i $pga_options $j");
}

sub build_three_consecutive_link_lirary{
    foreach(@three_consecutive_link){
        my $line = $_;
        my $times = grep{/$line/}@three_consecutive_link;
        if ( $times >= $min_coverage){
            push @three_consecutive_link_lirary,$line if ((grep{/$line/}@three_consecutive_link_lirary)==0);
        }
    }
}

sub connection_paths_evaluate{
    open( PgaResult, "$tmpdir/ConnectionPathSets.txt") || die "[Error]\n\tCannot find the file: ./$tmpdir/ConnectionPathSets.txt\n";
    my (@paths_library);
    my %Id_contigs = reverse %contigs_id;

    while (<PgaResult>) {
        next if !/^\d+/;
        chomp;
        my @forward_link = split /->/;
        my @reverse_link = reverse @forward_link;
        my ( $linkFr, $linkRv);
        for ( my $n = 1 ; $n <= 2 ; $n++ ) {
            for ( my $i = 0 ; $i < $#forward_link ; $i++ ) {
                $forward_link[$i] = sprintf( "%05d", $forward_link[$i] );
                $linkFr = $linkFr . $forward_link[$i]. "->";
            }
            for ( my $i = 1 ; $i <= $#reverse_link ; $i++ ) {
                $reverse_link[$i] = sprintf( "%05d", $reverse_link[$i] );
                $linkRv = $linkRv . $reverse_link[$i]. "->";
            }
        }
        push @paths_library,$linkFr;
        push @paths_library,$linkRv;
    }
    close PgaResult;

    my $maxscore = my $j = my $k = 0;
    my @contigs_link_score;
    foreach(@paths_library){
        my ( $sum_score, $link_score );
        my @contigs = split /->/;
        for ( my $i = 0 ; $i < $#contigs ; $i++ ) {
            my $ctgs_link = $contigs[$i] . "->" . $contigs[$i + 1];
            my $link_repeat_times = grep { /$ctgs_link/ } @paths_library;
            $link_score = $link_score . $link_repeat_times . "->";
            $sum_score += $link_repeat_times;
            my $orientation_detect = grep{/$Id_contigs{$contigs[$i]}/}grep{/$Id_contigs{$contigs[$i+1]}/} @orientation_library;
            $sum_score -= 10 if ($orientation_detect ==0 ); 
        }
        $link_score .= 0;
        push @contigs_link_score ,$link_score ;

        if ( $sum_score > $maxscore ) {
            $maxscore = $sum_score;
            $j        = $k;
        }
        $k++;
    }

    open(Repeat,">>$tmpdir/Repeat.txt");
    print Repeat "Possible repeat connection:(mincoverage:$min_coverage)\n";
    my ($path_link,$quality_score)=&restore_repeat_contigs($paths_library[$j],$contigs_link_score[$j]);
    close Repeat;
    open(OptimumPath,">>$tmpdir/OptimumPath.txt");
    print OptimumPath "Optimum contigs_Id connection path:\n$path_link\nOptimum contigs_Id connection score:\n$quality_score\n";
    close OptimumPath;

    my @link  = split /->/, $path_link;
    my @score = split /->/, $quality_score;
    my $num   = 1;
    my $min_score = $candidate_paths_num * $link_repetition_rate;
    my ($link_paths,$last,@paths,@scaffolded_paths);
    for ( my $i = 0 ; $i <= $#score ; $i++ ) {
        my ($orientation) = grep {/$Id_contigs{$link[$i+1]}/}grep {/^$Id_contigs{$link[$i]}/} @orientation_library;
        my ( $forward_orient, $backward_orient ) = $orientation =~ /\S+\s+\S+\s+(\S+)_(\S+)/;

        if ($forward_orient ne "F" && $forward_orient ne "R" ){
            $forward_orient = 1;
            $backward_orient = 0;
        }

        if ( $i == 0 && $score[0] >= $min_score ) {
            $link_paths = $link[0] . "_<$forward_orient>--->";
            $last  = $backward_orient;
            $num++;
            next;
        }

        if ( !defined($last) ) {
            $last = $forward_orient ;
            if ($forward_orient eq "F" || $forward_orient eq "R"){
                $last = $forward_orient ;
            }else{
                $last = $backward_orient;
            }
        }

        if ( $score[$i] >= $min_score && $last eq $forward_orient ) {
            $link_paths = $link_paths . $link[$i] . "_<$forward_orient>--->";
            $num++;
            $last = $backward_orient;
        }else {
            $link_paths = $num . "\t" . $link_paths . $link[$i] . "_<$last>";
            push @paths, $link_paths if (!( grep {/$link_paths/} @paths ));
            $link_paths = ();
            $num   = 1;
            $last = ();
        }
    }

    if ( @paths == 1 ) {
        my ($link) = $paths[0] =~ /\d+\s+(\S+)/;
        my @array = split /--->/, $link;
        my $contigs_link;
        for ( my $i = 0 ; $i < $match_contigs_num-1 ; $i++ ) {
            $contigs_link = $contigs_link . $array[$i] . "--->";
        }
        $contigs_link = $contigs_link . $array[-1];
        push @scaffolded_paths , $contigs_link;
    }else {
        foreach (@paths){
            my ($link) = $_ =~ /\d+\s+(\S+)/;
            my @array = grep {/$link/ } @paths;
            my @sorted_array = sort { $b <=> $a } @array;
            push @scaffolded_paths,$link if ( $_ eq $sorted_array[0] );
        }
    }

    open( Scaffold, ">$tmpdir/ScaffoldsPath.txt" ) || die;
    my $count = 1;
    foreach (@scaffolded_paths){
        my @contigs_pga = split /--->/;
        my $scaffold;
        for ( my $i = 0 ; $i <= $#contigs_pga ; $i++ ) {
            my ( $pga_num, $ctg_orient ) = $contigs_pga[$i] =~ /(\d+)_(<\S+>)$/;
            $scaffold = $scaffold . $Id_contigs{$pga_num} . $ctg_orient . "--->";
        }
        $scaffold=~s/--->$//;
        $count = sprintf("%06d",$count);
        print Scaffold ">scaffold_$count\n$scaffold\n";
        $count++;
    }
    close Scaffold;
    $totalNum_scaffolds = $count - 1 + $unmatch_contigs_num;
}

sub restore_repeat_contigs{
    my ($paths,$scores)=@_;
    my @path=split /->/,$paths;
    my @score=split /->/,$scores;
    my ($path_link,$quality_score);
    for(my $i=0;$i<=$#path;$i++){
        $path_link = $path_link.$path[$i]."->";
        $quality_score = $quality_score.$score[$i]."->";
        foreach(@three_consecutive_link_lirary){
            my ($left,$center,$right)=split /->/,$_;
            if($path[$i]==$contigs_id{$left} && $path[$i+1]==$contigs_id{$right}){
                print Repeat "\t$path[$i]"."->[$contigs_id{$center}]"."->$path[$i+1]\n";
                $path_link = $path_link.$contigs_id{$center}."->";
                $quality_score = $quality_score."10->";
            }
        }
    }
    return($path_link,$quality_score);
}

sub get_gap_fasta{
    open(Fastq,"$pacbio_reads");
    my ($i,$id,$start,$end);
    my $lable=0;
    while(<Fastq>){
        if($i % 4 == 0){
            my ($name)=$_=~/^@(\S+)/;
            if(my ($gap_information)=grep{/$name\t/}@gap_info){
                ($id,$start,$end)=$gap_information=~/^(\S+)\s+(\S+)\s+(\S+)\s+/;
                $lable=1;
            }
        }elsif($i%4 == 1 && $lable==1){
            my $len=$end-$start+1;
            my $gap_sequence=substr($_,$start,$len);
            $gap_fasta{$id}=$gap_sequence;
            $lable=0;
        }
        $i++;
    }
    close Fastq;
}

sub congtigs_scaffolding {
    open( Scaffold, "$tmpdir/ScaffoldsPath.txt" ) || die "[Error]\n\tCannot find the file: $tmpdir/ScaffoldsPath.txt\n";
    open( Gap, ">$tmpdir/Gap.fasta" );
    my $name;
    my $count = 1;
    while (<Scaffold>) {
        next if (/^>/);
        chomp;
        my @contigs  = split /--->/;
        my $numofcontigs = $#contigs + 1;
        if($numofcontigs == 1){
            ($name) = $_ =~ /^(\S+)<\S+>$/;
            $name = ">" . $name;
        }else{
            $count = sprintf("%06d",$count);
            $name = ">Scaffold_$count(Consists_of_$numofcontigs"."_contigs)";
            $count++;
        }
        my $scaffold = "$name\n";
        for ( my $n = 0 ; $n <= $#contigs ; $n++ ) {
            my ( $sequence, $gap_seq );
            my ( $contig, $orient ) = $contigs[$n] =~ /^(\S+)<(\S+)>$/;
            my $sequence = $contigs_seq{$contig};
            if ( $orient eq "R" ) {
                $sequence = reverse $sequence;
                $sequence =~ tr{wsatugcyrkmbdhvnATUGCYRKMBDHV\n}{WSTAACGRYMKVHDBNTAACGRYMKVHDB}d;
            }
            if ( $n == 0 ) {
                $scaffold .= $sequence;
                next;
            }
            my ($contig_previous) = $contigs[ $n - 1 ] =~ /^(\S+)<(\S+)>$/;
            my ($gap_len) = map {$_->[1] }map {[ $_, ( split /\s+/ )[2] ] } grep { /$contig/ } grep { /$contig_previous/ } @distance_library;
            print Gap ">$contig_previous\_$contig\n";

            if ( $gap_len >= 0 ) {
                if($gap_on_off eq "T"){
                    my @gap_array = grep{/$contig_previous/}grep{/$contig/}@gap_info;
                    my @sorted_gap_array=map {$_->[0] }sort {$b->[1] <=> $a->[1] } map {[ $_, ( split /\s+/ )[5] ] } @gap_array;
                    my $clustalo_coverage = @sorted_gap_array;
                    $clustalo_coverage=$max_gapfill_Reads if ($clustalo_coverage > $max_gapfill_Reads);
                    if($clustalo_coverage >= $min_gapfill_coverage){
                        open (TEMP,">$tmpdir/temp.fa");
                        for(my $j=0;$j<$clustalo_coverage;$j++){
                            my ($fastq_name)=$sorted_gap_array[$j]=~/^(\S+)\s+/;
                            if($sorted_gap_array[$j]=~/$contig_previous\t$contig/){
                                print TEMP ">$j\n$gap_fasta{$fastq_name}\n";
                            }elsif($sorted_gap_array[$j]=~/$contig\t$contig_previous/){
                                my $sequence_gap = reverse $gap_fasta{$fastq_name};
                                $sequence_gap =~ tr{wsatugcyrkmbdhvnATUGCYRKMBDHV\n}{WSTAACGRYMKVHDBNTAACGRYMKVHDB}d;
                                print TEMP ">$j\n$sequence_gap\n";
                            }
                        }
                        close TEMP;
                        system("clustalo -i $tmpdir/temp.fa -o $tmpdir/clustalo.fa -t DNA --threads $threads --force");
                        $gap_seq = &gap_filling;
                        unlink "$tmpdir/temp.fa";
                        unlink "$tmpdir/clustalo.fa";
                    }else{
                        $gap_seq = "N" x $gap_len;
                    }
                }else{
                    $gap_seq = "N" x $gap_len;
                }
                print Gap "$gap_seq\n";
            }else {
                $gap_len  = abs($gap_len) + 3;
                $sequence = substr( $sequence, $gap_len );
                $gap_seq  = "NNN";
                print Gap "$gap_seq\t(Overlapped)\n";
            }
            $scaffold = $scaffold . $gap_seq . $sequence;
        }
        $scaffold .= "\n";
        my $string = IO::String->new($scaffold);
        my $seqio = Bio::SeqIO->new( -fh => $string, -format => 'fasta' );
        my $out = Bio::SeqIO->new( -file => ">>$outputfile", -format => 'fasta' );
        while ( my $seq = $seqio->next_seq ) {
            $out->write_seq($seq);
        }
    }
    if(-e "$tmpdir/PacBio.unmatched.fasta"){
        system("mv $outputfile $tmpdir/scaffolds.temp");
        system("cat $tmpdir/scaffolds.temp $tmpdir/PacBio.unmatched.fasta > $outputfile");
        unlink "$tmpdir/scaffolds.temp";
    }
    close Gap;
    close Scaffold;
}

sub gap_filling{
    my $inClutalwOut = Bio::SeqIO->new(-file => "$tmpdir/clustalo.fa", -format => 'fasta');
    my ($row,$column);
    my $j=0;
    my @AoA;
    while ( my $seq = $inClutalwOut->next_seq() ) {
        my $seqence = $seq->seq;
        $seqence=~s/A/1 /g;
        $seqence=~s/G/2 /g;
        $seqence=~s/C/3 /g;
        $seqence=~s/T/4 /g;
        $seqence=~s/-/0 /g;
        my @a=split / /,$seqence;
        $column = $#a;
        for(my $i=0;$i<@a;$i++){
            $AoA[$i][$j]=$a[$i];
        }
        $row=$j;
        $j++;
    }

    my @c;
    for(my$i=0;$i<=$column;$i++){
        my @b;
        for(my$j=0;$j<=$row;$j++){
            push @b,$AoA[$i][$j];
        }
        my $stat = Statistics::Descriptive::Full->new();
        $stat->add_data( \@b );
        my $mode = $stat->mode();
        push @c,$mode;
    }

    my $string = join('',@c);
    $string=~s/1/a/g;
    $string=~s/2/g/g;
    $string=~s/3/c/g;
    $string=~s/4/t/g;
    $string=~s/0//g;

    return $string;
}

sub statistics{
    my($file)=@_;
    my($Total_length,$length,$Largest_contig,$N50,$N90,$A,$G,$C,$T);
    my (@x,@statistics);
    open(File,"$file");
    while(<File>){
        if(/^[\>\@]/){
            if($length>0){
                $Total_length += $length;
                push@x,$length;
            }
            $length = 0;
        }else{
            s/\s//g;
            $length += length($_);
            my $num_a = $_ =~ tr/Aa//;
            $A += $num_a;
            my $num_g = $_ =~ tr/Gg//;
            $G += $num_g;
            my $num_c = $_ =~ tr/Cc//;
            $C += $num_c;
            my $num_t = $_ =~ tr/Tt//;
            $T += $num_t;
            $Largest_contig = $length if ($length > $Largest_contig);
        }
    }
    if ($length > 0){
        $Total_length += $length;
        push @x,$length;
    }

    push @statistics,sprintf("%12d",$Total_length);
    push @statistics,sprintf("%12d",$Largest_contig);
    my $Average_length = sprintf("%12d",$Total_length / @x);
    push @statistics,$Average_length;
    my $GC = sprintf( "%12.2f", ($G + $C) / ($A + $G + $C + $T) * 100 );
    push @statistics,$GC;

    if(@x==1){
        # N50
        push @statistics,sprintf("%12d",$Total_length);
        # N90
        # push @statistics,sprintf("%12d",$Total_length);
    }else{
        @x=sort{$b<=>$a}@x;
        my ($count,$half);
        for (my $j=0;$j<@x;$j++){
            $count += $x[$j];
            if(($count >= $Total_length/2)&&($half==0)){
                # N50
                push @statistics,sprintf("%12d",$x[$j]);
                $half = $x[$j];
            #}elsif($count >= $Total_length * 0.9){
                # N90
                # push @statistics,sprintf("%12d",$x[$j]);
                # last;
            }
        }
    }
    close File;
    return (@statistics);
}

