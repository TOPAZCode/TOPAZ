#!/bin/bash


Here=$PWD

WorkDir=./LHC13_53_2.18



# ---------------------------------------------
# various grids for couplings

coupllist0=(+0.00)

coupllist1=(+0.00
           -0.50
           -0.10
           +0.10
           +0.50
          )

coupllist2=(+0.00
           -0.50
           -0.40
           -0.30
           -0.20
           -0.10
           +0.10
           +0.20
           +0.30
           +0.40
           +0.50
          )

coupllist3=(+0.00
           -0.50
           -0.45
           -0.40
           -0.35
           -0.30
           -0.25
           -0.20
           -0.15
           -0.10
           -0.05
           +0.05
           +0.10
           +0.15
           +0.20
           +0.25
           +0.30
           +0.35
           +0.40
           +0.45
           +0.50
          )
# ---------------------------------------------
# select two of the above grids

coupllistV=( ${coupllist2[*]} );

coupllistA=( ${coupllist3[*]} );



# ---------------------------------------------

cd $WorkDir
mkdir ./Results
touch ./Results/Summary.dat
echo "#Summary file" >  ./Results/Summary.dat
for couplV in "${coupllistV[@]}"
do
for couplA in "${coupllistA[@]}"
do

  summer addsilent ./71/LHC.71.LOLO_V${couplV}A${couplA}.dat ./72/LHC.72.LOLO_V${couplV}A${couplA}.dat ./Results/LHC.LOLO_V${couplV}A${couplA}.dat
  
  TotalCS=`summer sumbins ./Results/LHC.LOLO_V${couplV}A${couplA}.dat 16 | tail -1`
  echo ${couplV} ${couplA} $TotalCS 
  echo ${couplV} ${couplA} $TotalCS >> ./Results/Summary.dat

  gnuplotStr+="'LHC.LOLO_V${couplV}A${couplA}.dat' using (\$2):( \$1==NHi ? \$3*8*300/${TotalCS} : 1/0) notitle w st, "
#   read
done
done


# ---------------------------------------------
# this is for printing the gnuplot string
#  echo $gnuplotStr







