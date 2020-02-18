#/bin/bash

if [[ "$1" =~ ".seq" ]]; then 
  echo ------ doing: $0 $1
else
  echo do: $0 _FILE_NAME_.seq -e \(OR -Emissivity\) 0.95
  break
fi

setEmissivity=
if [[ "$2" = "-e" || "$2" = "-Emissivity"  ]] && [[ "$3" != "" ]]; then 
  setEmissivity=$3
fi


binarydir=fout
rm -rf $binarydir
mkdir -p $binarydir

echo '------ converting seq to binary ------ ' 
#
# read the sequence file
# produce the binary files frame by frame as seq_n*.fff
#

./share/seqToBin.py $1 $binarydir     #This should work for any number of sequence files...
#./share/seqtobinary.pl $1 $binarydir #This method works for small sequence files
echo ''

txtoutdir=tout
rm -rf $txtoutdir
mkdir -p $txtoutdir
echo '------ converting binary to text ----- '
./share/binarytotext.sh $binarydir $txtoutdir $setEmissivity
echo ''

nfile=`ls -l $binarydir/*_*.* | wc -l`
echo '------ converting text to root ----- '
echo '    found number of files: '$nfile''
outdir=roo
rm -rf $outdir
mkdir -p $outdir
if (( nfile <= 0 )); then
  echo 'number of files '$nfile' <= 0 '
  break
fi
./share/texttoroot.py $outdir $nfile
echo ''


echo 'clear binary folder: '$binarydir''
rm -rf $binarydir
echo 'clear text out folder: '$txtoutdir''
rm -rf $txtoutdir

mv config $outdir
ls -l $outdir/config

echo 'root results in folder: '$outdir''
echo ''
echo 'All done!'

