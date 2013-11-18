#! /bin/bash
#
#
#
CollTag="LHC13"
#
#
#
#
echo "runTOPAZ: Export general setup"
. /afs/cern.ch/sw/lcg/hepsoft/0.8/x86_64-slc6-gcc47-opt/setup.sh
export myAFSpath=/afs/cern.ch/user/m/maschulz/batch/TOPAZ


echo "runTOPAZ: Copying executable and pdfs to cluster"
cp ${myAFSpath}/TOPAZ .
mkdir ./PDFS/
cp ${myAFSpath}/PDFS/* ./PDFS/


# create directory at home
echo "runTOPAZ: Create directories at home for $CollTag"
mkdir ${myAFSpath}/${CollTag}_${ObsSet}_${Mu}
mkdir ${myAFSpath}/${CollTag}_${ObsSet}_${Mu}/${Process}




# convert wall time from HH:M into seconds
IFS=: read h m <<< "$WALLTIME"
WT_SEC=$(( 10#$h * 3600 + 10#$m * 60 ))
# define half walltime and close-to-walltime (=180 seconds before)
let "WT_SEC_HALF= $WT_SEC /2"
let "WT_SEC_ALMOST=$WT_SEC - 180"


echo "runTOPAZ: full   walltime = $WT_SEC sec"
echo "runTOPAZ: half   walltime = $WT_SEC_HALF sec"
echo "runTOPAZ: almost walltime = $WT_SEC_ALMOST sec"


# running the executable
./TOPAZ ObsSet=$ObsSet Collider=$Collider Process=$Process Correction=$Correction NLOParam=$NLOParam PDFSet=$PDFSet ZDK=$ZDK TopDK=$TopDK VegasIt0=$VegasIt0 VegasNc0=$VegasNc0 VegasIt1=$VegasIt1 VegasNc1=$VegasNc1 MuRen=$Mu MuFac=$Mu FileTag=$FileTag VegasSeed=$VegasSeed DipAlpha=$DipAlpha DipAlpha2=$DipAlpha2 DKAlpha=$DKAlpha GridIO=$GridIO HelSamp=$HelSamp GridFile=$GridFile RelDelF1A=$RelDelF1A RelDelF1V=$RelDelF1V &

TOPAZPID=$!
STARTTIME=$(date +%s)
echo `date` "(" ${STARTTIME} ")"
echo "runTOPAZ: Execute TOPAZ in background (PID=${TOPAZPID})"


while ps -p $TOPAZPID >/dev/null
do
#   check every 2 min. if TOPAZ is still running
    sleep 60
    CURRTIME=$(date +%s)
    RUNTIME=$((CURRTIME-STARTTIME))
    echo "runTOPAZ: Runtime check ${RUNTIME}"
    if [ $RUNTIME -gt $WT_SEC_HALF ] &&  [ $RUNTIME -lt $((WT_SEC_HALF+80)) ]; then
        cp ./${CollTag}_${ObsSet}_${Mu}/${Process}/* ${myAFSpath}/${CollTag}_${ObsSet}_${Mu}/${Process}
        date
        echo "runTOPAZ: Copied preliminary results back home"
    fi 
    if [ $RUNTIME -gt $WT_SEC_ALMOST ] &&  [ $RUNTIME -lt $((WT_SEC_ALMOST+80)) ]; then
        cp ./${CollTag}_${ObsSet}_${Mu}/${Process}/* ${myAFSpath}/${CollTag}_${ObsSet}_${Mu}/${Process}
        date
        echo "runTOPAZ: Copied preliminary results back home"
    fi 
done


# waiting for TOPAZ (in background) to finish
date
CURRTIME=$(date +%s)
RUNTIME=$((CURRTIME-STARTTIME))
echo "runTOPAZ: TOPAZ finished after ${RUNTIME} seconds"


# copying final files from cluster to home
#
cp ./${CollTag}_${ObsSet}_${Mu}/${Process}/* ${myAFSpath}/${CollTag}_${ObsSet}_${Mu}/${Process}
echo "runTOPAZ: Copied final results back home"

