#!/bin/bash                                                          
clear                                                                


Multiplier="*(8.0)";


WorkDir=./LHC13_53_2.18_Test02/

RelDelF1A=+0.00
RelDelF1V=+0.00
FileTag="_V${RelDelF1V}A${RelDelF1A}"

                                             
filelistLOLO=( 71/LHC.71.LOLO${FileTag}.dat
               72/LHC.72.LOLO${FileTag}.dat
              )

filelistNLO1=( 71/LHC.71.LO.dat
               72/LHC.72.LO.dat
               71/LHC.71.1L_070.dat
               72/LHC.72.1L_070.dat
               71/LHC.71.DK1L0.dat
               71/LHC.71.DKRE0.dat
               72/LHC.72.DK1L0.dat
               72/LHC.72.DKRE0.dat
               73/LHC.73.ID1110.dat
               73/LHC.73.RE1110.dat
               77/LHC.77.ID1110.dat
               77/LHC.77.RE1110.dat
               75/LHC.75.ID1110.dat
               75/LHC.75.RE1110.dat
               76/LHC.76.ID1110.dat
               76/LHC.76.RE1110.dat)                            

filelistNLO2=( 71/LHC.71.LO.dat
               72/LHC.72.LO.dat
               71/LHC.71.1L_070.dat
               72/LHC.72.1L_070.dat
               71/LHC.71.DK1L1.dat
               71/LHC.71.DKRE1.dat
               72/LHC.72.DK1L1.dat
               72/LHC.72.DKRE1.dat
               73/LHC.73.ID2220.dat
               73/LHC.73.RE2220.dat
               77/LHC.77.ID2220.dat
               77/LHC.77.RE2220.dat
               75/LHC.75.ID2220.dat
               75/LHC.75.RE2220.dat
               76/LHC.76.ID2220.dat
               76/LHC.76.RE2220.dat)                    



filelist=( "${filelistLOLO[@]}" )
Counter=0
Sum="0"
SumError="0"
for file in "${filelist[@]}"
do
    let Counter=$Counter+1
    TotCS[$Counter]=-999999999
    Error[$Counter]=-999999999
    TotCS[$Counter]=`grep -e "TotCS" ./$WorkDir/$file | cut -c12-35 | awk -F"E" 'BEGIN{OFMT="%30.16f"} {print $1 * (10 ^ $2)}'`
    Error[$Counter]=`grep -e "TotCS" ./$WorkDir/$file | cut -c36-60 | awk -F"E" 'BEGIN{OFMT="%30.16f"} {print $1 * (10 ^ $2)}'`
     echo ""
     echo "file($Counter): $file"
     echo "TotCS: ${TotCS[$Counter]}"
     echo "Error: ${Error[$Counter]}"
    Sum=$Sum" + (${TotCS[$Counter]})  ${Multiplier}"
    SumError=$SumError" + (${Error[$Counter]} ${Multiplier} )^2"
done
echo
echo "-------------------------------------------"
echo ""
echo "LOLO"
 echo "Sum = " $Sum
export Sum=`echo $Sum | bc -l`
# printf "Sum = $Sum \n"
printf 'Sum      = %16.8e \n' $Sum
# echo "SumError = sqrt[" $SumError "]"
export SumError=`echo "sqrt($SumError)" | bc -l`
# printf "SumError = $SumError \n"
printf 'SumError = %16.8e \n' $SumError
echo ""

read


filelist=( "${filelistNLO1[@]}" )
Counter=0
Sum="0"
SumError="0"
for file in "${filelist[@]}"
do
    let Counter=$Counter+1
    TotCS[$Counter]=-999999999
    Error[$Counter]=-999999999
    TotCS[$Counter]=`grep -e "TotCS" ./$WorkDir/$file | cut -c12-35 | awk -F"E" 'BEGIN{OFMT="%30.16f"} {print $1 * (10 ^ $2)}'`
    Error[$Counter]=`grep -e "TotCS" ./$WorkDir/$file | cut -c36-60 | awk -F"E" 'BEGIN{OFMT="%30.16f"} {print $1 * (10 ^ $2)}'`
#     echo ""
#     echo "file($Counter): $file"
#     echo "TotCS: ${TotCS[$Counter]}"
#     echo "Error: ${Error[$Counter]}"
    Sum=$Sum" + (${TotCS[$Counter]})"
    SumError=$SumError" + (${Error[$Counter]})^2"
done
echo "NLO1"
# echo "Sum = " $Sum
export Sum=`echo $Sum | bc -l`
# printf "Sum = $Sum \n"
printf 'Sum      = %16.8e \n' $Sum
# echo "SumError = sqrt[" $SumError "]"
export SumError=`echo "sqrt($SumError)" | bc -l`
# printf "SumError = $SumError \n"
printf 'SumError = %16.8e \n' $SumError
echo ""




filelist=( "${filelistNLO2[@]}" )
Counter=0
Sum="0"
SumError="0"
for file in "${filelist[@]}"
do
    let Counter=$Counter+1
    TotCS[$Counter]=-999999999
    Error[$Counter]=-999999999
    TotCS[$Counter]=`grep -e "TotCS" ./$WorkDir/$file | cut -c12-35 | awk -F"E" 'BEGIN{OFMT="%30.16f"} {print $1 * (10 ^ $2)}'`
    Error[$Counter]=`grep -e "TotCS" ./$WorkDir/$file | cut -c36-60 | awk -F"E" 'BEGIN{OFMT="%30.16f"} {print $1 * (10 ^ $2)}'`
#     echo ""
#     echo "file($Counter): $file"
#     echo "TotCS: ${TotCS[$Counter]}"
#     echo "Error: ${Error[$Counter]}"
    Sum=$Sum" + (${TotCS[$Counter]})"
    SumError=$SumError" + (${Error[$Counter]})^2"
done
echo "NLO2"
# echo "Sum = " $Sum
export Sum=`echo $Sum | bc -l`
# printf "Sum = $Sum \n"
printf 'Sum      = %16.8e \n' $Sum
# echo "SumError = sqrt[" $SumError "]"
export SumError=`echo "sqrt($SumError)" | bc -l`
# printf "SumError = $SumError \n"
printf 'SumError = %16.8e \n' $SumError
echo 


