#!/usr/bin/perl

use strict;
use warnings;
use utf8;
use Cwd "abs_path";

# centre of mass for each of the seeds

## these seed coordinates are in LPI order
my %seedCOMs = (
		"R_PPC"        =>  [ 44, -50,  50],
		"R_PPC_SPL"    =>  [ 46, -54,  46], 
		"L_PPC_SPL"    =>  [-42, -52,  48], 
		"R_hipp_MTG"   =>  [ 38, -30, -06], ## white matter 4mm from hippocampus
		"MPFC_1"       =>  [-02,  38,  12], ## perigenual cingulate
		"L_DLPFC_1"    =>  [-42,  26,  32], 
		"L_DLPFC_2"    =>  [-38,  26,  36], 
		"L_DLPFC_3"    =>  [-34,  34,  36],
		"Mid_cing"     =>  [-02,  04,  22], ## midcingulate/corpus callosum
		"MPFC_2"       =>  [-02,  46,  16], ## perigenual cingulate
		"R_Precuneus"  =>  [ 18, -66,  34]
	       );

my $seedFile=abs_path("../data/config/kasier_supplemental_seeds.txt");
if ( -f $seedFile ) {
  unlink $seedFile;
}

print "Opening $seedFile\n";
open (my $SH, '>>', $seedFile) || die "Unable to open $seedFile\n";

foreach my $seed ( sort keys %seedCOMs ) {
  my @com=@{$seedCOMs{$seed}};
  # 3dUndump expects a file of xyz coordinates, so tell it to read from standard input and echo the coordinates instead.
  my $command="echo @com | 3dUndump -master ../standard/MNI152_T1_3mm_brain.nii.gz -srad 3 -orient LPI -prefix " . abs_path("../data/seeds/kaiser") . "/$seed -xyz -";
  print "$command\n";
  system $command;

  my $seedFileLine="\$DATA/seeds/kaiser/$seed+tlrc.HEAD\n";
  print $SH $seedFileLine;
  #$command="3dNotes -a \"Coordinates are @com\" " . abs_path("../data/restingmasks") . "/acc$seed.nii.gz";
  #print "$command\n";
  #system ($command);
}

close $SH;
