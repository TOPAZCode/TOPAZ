#!/bin/bash


RunWhat=$1

export RunDir=$PWD

export ObsSet=56
export TopDK=4
export ZDK=1
#
export Collider=13
export PDFSet=1
#
#export Mu=1.09
export Mu=2.18
#export Mu=4.36
#
#
#
export VegasSeed=19
export DipAlpha=00000 DipAlpha2=1 DKAlpha=111
export FileTag="."
export GridFile="."
export HelSamp=0 GridIO=0
#
#
#
# ---------------------------------------------
coupllist0=(+0.00 )
coupllist1=(+0.00-0.50 -0.10 +0.10 +0.50 )
coupllist2=(+0.00 -0.50 -0.40 -0.30 -0.20 -0.10 +0.10 +0.20 +0.30 +0.40 +0.50 )
coupllist2a=(-1.00 -0.90 -0.80 -0.70 -0.60 +0.60 +0.70 +0.80 +0.90 +1.00      )
coupllist3=(+0.00 -0.50 -0.45 -0.40 -0.35 -0.30 -0.25 -0.20 -0.15 -0.10 -0.05 
            +0.05 +0.10 +0.15 +0.20 +0.25 +0.30 +0.35 +0.40 +0.45 +0.50       )
coupllist3a=(+0.025 -0.475 -0.425 -0.375 -0.325 -0.275 -0.225 -0.175 -0.125 -0.075 
             -0.025 +0.075 +0.125 +0.175 +0.225 +0.275 +0.325 +0.375 +0.425 +0.475 )
coupllist4=(-0.10)
coupllist5=(-0.10 -0.05 +0.05 +0.10 +0.15 )



coupllistV1=(-1.40 -1.20 -1.00 -0.80 -0.60 -0.40 -0.20 +0.00 +0.20 +0.40 +0.60 +0.80 +1.00 +1.20 +1.40) 
coupllistA1=(-0.20 -0.15 -0.10 -0.05 +0.00 +0.05 +0.10 +0.15 +0.20)

coupllistV1scale=(-1.00 +1.00) 
coupllistA1scale=(-0.10 +0.10)


# ---------------------------------------------
coupllistV=( ${coupllist0[*]} );
coupllistA=( ${coupllist0[*]} );

echo Ready to run the job: $RunWhat ?
read


# ---------------------------------------------


for couplV in "${coupllistV[@]}"
do
for couplA in "${coupllistA[@]}"
do

#if [ "$couplV" == "+0.00" ] && [ "$couplA" == "+0.00" ] ; then
#   echo "Skip the (+0.00,+0.00) contribution"
#   continue
#fi

if [ "$RunWhat" == "LOLO" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=02:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export FileTag="_V${RelDelF1V}A${RelDelF1A}"
    export Correction=0 NLOParam=1
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=1000000 VegasNc1=2000000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=1000000 VegasNc1=2000000
    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
fi




# ---------------------------------------------


if [ "$RunWhat" == "LO" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=02:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export FileTag="_V${RelDelF1V}A${RelDelF1A}"
    export Correction=0 NLOParam=2
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=1000000 VegasNc1=2000000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=1000000 VegasNc1=2000000
    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
fi


# ---------------------------------------------


if [ "$RunWhat" == "VI" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=48:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export FileTag="_V${RelDelF1V}A${RelDelF1A}"
    export Correction=1 NLOParam=2 GridIO=2
    export VegasSeed=20
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

    export VegasSeed=40
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

    export VegasSeed=60
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

    export VegasSeed=80
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=30000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

    export VegasSeed=19
fi




# ---------------------------------------------




if [ "$RunWhat" == "RE" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=48:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export Correction=2 NLOParam=2
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
#    export DipAlpha=22220 DipAlpha2=5 FileTag="1111_V${RelDelF1V}A${RelDelF1A}"
#    export Process=73 VegasIt0=3 VegasIt1=5 VegasNc0=10000000 VegasNc1=10000000
#    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=74 VegasIt0=3 VegasIt1=5 VegasNc0=10000000 VegasNc1=10000000
#    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=75 VegasIt0=3 VegasIt1=5 VegasNc0=10000000 VegasNc1=10000000
#    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=76 VegasIt0=3 VegasIt1=5 VegasNc0=10000000 VegasNc1=10000000
#    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh


    export DipAlpha=22220 DipAlpha2=1 FileTag="2222_V${RelDelF1V}A${RelDelF1A}"
    export Process=73 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=4000000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=74 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=4000000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=75 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=4000000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=76 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=4000000
    bsub -q 2nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
fi



# ---------------------------------------------




if [ "$RunWhat" == "ID" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=24:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export Correction=3 NLOParam=2
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
#    export DipAlpha=22220 DipAlpha2=5 FileTag="1111_V${RelDelF1V}A${RelDelF1A}"
#    export Process=73 VegasIt0=3 VegasIt1=5 VegasNc0=5000000 VegasNc1=5000000
#    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=74 VegasIt0=3 VegasIt1=5 VegasNc0=5000000 VegasNc1=5000000
#    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=75 VegasIt0=3 VegasIt1=5 VegasNc0=5000000 VegasNc1=5000000
#    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=76 VegasIt0=3 VegasIt1=5 VegasNc0=5000000 VegasNc1=5000000
#    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh


    export DipAlpha=22220 DipAlpha2=1 FileTag="2222_V${RelDelF1V}A${RelDelF1A}"
    export Process=73 VegasIt0=3 VegasIt1=5 VegasNc0=2000000 VegasNc1=2000000
    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=74 VegasIt0=3 VegasIt1=5 VegasNc0=2000000 VegasNc1=2000000
    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=75 VegasIt0=3 VegasIt1=5 VegasNc0=2000000 VegasNc1=2000000
    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=76 VegasIt0=3 VegasIt1=5 VegasNc0=2000000 VegasNc1=2000000
    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
fi




# ---------------------------------------------




if [ "$RunWhat" == "DKVI" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=08:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export Correction=4 NLOParam=2
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    export DipAlpha=00001 DKAlpha=111 FileTag="1_V${RelDelF1V}A${RelDelF1A}"
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=2000000 VegasNc1=2000000
    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=2000000 VegasNc1=2000000
    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh


#    export DipAlpha=00001 DKAlpha=555 FileTag="5_V${RelDelF1V}A${RelDelF1A}"
#    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=5000000 VegasNc1=5000000
#    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=5000000 VegasNc1=5000000
#    bsub -q 8nh -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
fi




# ---------------------------------------------




if [ "$RunWhat" == "DKRE" ] || [ "$RunWhat" == "ALL" ] ; then
    export WALLTIME=24:00
    export RelDelF1A=$couplA
    export RelDelF1V=$couplV
    export Correction=5 NLOParam=2
    export DipAlpha=00001 DKAlpha=111 export FileTag="1_V${RelDelF1V}A${RelDelF1A}"
    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=4000000 VegasNc1=4000000
    echo "Submit $RunWhat for couplings deltaV=$RelDelF1V deltaA=$RelDelF1A"
    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=400000 VegasNc1=4000000
    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh

#    export DipAlpha=00001 DKAlpha=555 export FileTag="5_V${RelDelF1V}A${RelDelF1A}"
#    export Process=71 VegasIt0=3 VegasIt1=5 VegasNc0=10000000 VegasNc1=10000000
#    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
#    export Process=72 VegasIt0=3 VegasIt1=5 VegasNc0=10000000 VegasNc1=10000000
#    bsub -q 1nd -J "${ObsSet}.${RunWhat}.${Process}" -n 1,1 -R "type==SLC6_64" -W $WALLTIME < runTOPAZ.sh
fi






if [ "$RunWhat" == "test" ] ; then
      echo "nothing to be done here"
fi



done
done

