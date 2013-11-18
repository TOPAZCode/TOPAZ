#! /bin/bash
#
#
#
CollTag="LHC13"
#
#
# general setup
. /afs/cern.ch/sw/lcg/hepsoft/0.8/x86_64-slc6-gcc47-opt/setup.sh
export myAFSpath=/afs/cern.ch/user/m/maschulz/batch/TOPAZ


echo "runTOPAZ: Copying executable and pdfs to cluster"
cp ${myAFSpath}/TOPAZ .
mkdir ./PDFS/
cp ${myAFSpath}/PDFS/* ./PDFS/



echo "runTOPAZ: Running TOPAZ. DATE = `date`"
./TOPAZ ObsSet=$ObsSet Collider=$Collider Process=$Process Correction=$Correction NLOParam=$NLOParam PDFSet=$PDFSet ZDK=$ZDK TopDK=$TopDK VegasIt0=$VegasIt0 VegasNc0=$VegasNc0 VegasIt1=$VegasIt1 VegasNc1=$VegasNc1 MuRen=$Mu MuFac=$Mu FileTag=$FileTag VegasSeed=$VegasSeed DipAlpha=$DipAlpha DipAlpha2=$DipAlpha2 DKAlpha=$DKAlpha GridIO=$GridIO HelSamp=$HelSamp GridFile=$GridFile RelDelF1A=$RelDelF1A RelDelF1V=$RelDelF1V DataDir=/afs/cern.ch/user/m/maschulz/batch/TOPAZ/

echo "runTOPAZ: Done. DATE = `date`"


