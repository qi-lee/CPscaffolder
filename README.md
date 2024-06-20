# CPscaffolder


### Basic usage : perl CPscaffolder.pl -i contigs -r PacBioReads [-options]

###  Example : 
       1. $ perl CPscaffolder.pl -i contigs.fasta -r pacbioreads.bas.h5
      2. $ perl CPscaffolder.pl -i contigs.fasta -b dir/blasr.m5 
       3. $ perl CPscaffolder.pl -i contigs.fasta -r pacbioreads.fastq -o scaffolds.fasta
       4. $ perl CPscaffolder.pl -i contigs.fasta -r pacbioreads.fastq -t 24 -g T 
       5. $ perl dir/CPscaffolder.pl -i dir/contigs.fa -r dir/pacbioreads.fq -c dir/config.ini 

###  Option Description (default_value) :
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

###   Software Requirement :
       - BLASR (http://bix.ucsd.edu/projects/blasr/)
       - BioPerl (http://www.bioperl.org/wiki/Installing_BioPerl)
       - Statistics::Descriptive Perl module (http://www.cpan.org/modules/INSTALL.html)
       - Config::IniFiles Perl module (http://www.cpan.org/modules/INSTALL.html)
