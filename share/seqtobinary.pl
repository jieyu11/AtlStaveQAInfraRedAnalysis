#!/usr/bin/perl

# @purpose:
#   Read .seq files from Flir IR camera and write each frame to temporary binary file.
#
# @usage:
#   seqtobinary.pl _FILE_NAME_.seq
#
# @note:
#   When first using this code for a new camera, it might need find the bits separating
#   each frame, which is possibly IR camera specific. Please run:
#     hexdump -n16 -C _FILE_NAME_.seq 
#
#   @@Example
#     >$ hexdump -n16 -C Rec-000667_test.seq 
#     00000000  46 46 46 00 52 65 73 65  61 72 63 68 49 52 00 00  |FFF.ResearchIR..|
#     00000010
#     So, for this camera, the separation patten is:
#     \x46\x46\x46\x00\x52\x65\x73\x65\x61\x72\x63\x68\x49\x52
#     which == FFFResearchIR
#

my $input = $ARGV[0];
if (not defined $input) {
  die "Please provide _FILE_NAME_.seq as argument! \n";
}

my $outfold = $ARGV[1];
if (not defined $outfold) {
  print "Output folder not provided. Using fout as default!\n";
  $outfold = fout;
}


undef $/;
$_ = <>;
$n = 0;

$pat="\x46\x46\x46\x00\x52\x65\x73\x65\x61\x72\x63\x68\x49\x52";

print "Input name: $input\n";
print "Frame separation patten: $pat\n";

for $content (split(/(?=$pat)/)) {
    open(OUT, ">$outfold/frame_" . $n . ".fff");
    binmode OUT;
    print OUT $content;
    close(OUT);
    ++$n;
}
