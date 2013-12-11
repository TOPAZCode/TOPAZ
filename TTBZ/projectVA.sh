#!/bin/bash

WorkDir=./LHC13_53_2.18


# ------------- LO data set ---------------------
# dV1=-1.40
# dA1=-0.20
# file1=$WorkDir/LHC.LOLO_V-1.40A-0.20.dat
# 
# dV2=-1.40
# dA2=+0.20
# file2=$WorkDir/LHC.LOLO_V-1.40A+0.20.dat
# 
# dV3=+1.40
# dA3=+0.20
# file3=$WorkDir/LHC.LOLO_V+1.40A+0.20.dat
# 
# dV4=+1.40
# dA4=-0.20
# file4=$WorkDir/LHC.LOLO_V+1.40A-0.20.dat
# 
# dV5=+0.00
# dA5=+0.00
# file5=$WorkDir/LHC.LOLO_V+0.00A+0.00.dat
# 
# dV6=+1.00
# dA6=+0.10
# file6=$WorkDir/LHC.LOLO_V+1.00A+0.10.dat


# ------------- NLO data set ---------------------
dV1=-1.40
dA1=-0.20
file1=$WorkDir/LHC.NLO2_V-1.40A-0.20.dat

dV2=-1.20
dA2=+0.15
file2=$WorkDir/LHC.NLO2_V-1.20A+0.15.dat

dV3=+1.40
dA3=+0.20
file3=$WorkDir/LHC.NLO2_V+1.40A+0.20.dat

dV4=+1.40
dA4=-0.10
file4=$WorkDir/LHC.NLO2_V+1.40A-0.10.dat

dV5=+0.00
dA5=+0.00
file5=$WorkDir/LHC.NLO2_V+0.00A+0.00.dat

dV6=+1.20
dA6=+0.10
file6=$WorkDir/LHC.NLO2_V+1.20A+0.10.dat


# ----------------------------------------------


NHi=17
outfiletag=NLO
summer ttbzcoupl $NHi $dV1 $dA1 $file1 $dV2 $dA2 $file2 $dV3 $dA3 $file3 $dV4 $dA4 $file4 $dV5 $dA5 $file5 $dV6 $dA6 $file6 $outfiletag




