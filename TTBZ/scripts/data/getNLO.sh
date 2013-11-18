#!/bin/bash
#
#
#
ThisDir=$pwd
WorkDir=./LHC13_56_4.36/
cd $WorkDir
#
#
RelDelF1A=+0.00
RelDelF1V=+0.00
FileTag="_V${RelDelF1V}A${RelDelF1A}"
#
#
summer add ./71/LHC.71.LOLO${FileTag}.dat ./72/LHC.72.LOLO${FileTag}.dat ./LHC.LOLO${FileTag}.dat
read
#
summer add ./71/LHC.71.LO${FileTag}.dat ./72/LHC.72.LO${FileTag}.dat ./LHC.LO${FileTag}.dat
read
#
#
#
#
#summer avg ./71/LHC.71.1L_020${FileTag}.dat  ./71/LHC.71.1L_070${FileTag}.dat  ./71/LHC.71.1L_060${FileTag}.dat ./71/LHC.71.1L_080${FileTag}.dat ./LHC.71.1L${FileTag}.dat
#read
#
#summer avg ./72/LHC.72.1L_020${FileTag}.dat  ./72/LHC.72.1L_070${FileTag}.dat  ./72/LHC.72.1L_060${FileTag}.dat ./72/LHC.72.1L_080${FileTag}.dat ./LHC.72.1L${FileTag}.dat
#read
#
summer add ./71/LHC.71.1L${FileTag}.dat ./72/LHC.72.1L${FileTag}.dat ./LHC.1L${FileTag}.dat
read
#
#
# alpha=low
summer add ./73/LHC.73.ID1111${FileTag}.dat ./73/LHC.73.RE1111${FileTag}.dat ./LHC.73.REID1111${FileTag}.dat
read
#
summer add ./74/LHC.74.ID1111${FileTag}.dat ./74/LHC.74.RE1111${FileTag}.dat ./LHC.74.REID1111${FileTag}.dat
read
#
summer add ./75/LHC.75.ID1111${FileTag}.dat ./75/LHC.75.RE1111${FileTag}.dat ./LHC.75.REID1111${FileTag}.dat
read
#
summer add ./76/LHC.76.ID1111${FileTag}.dat ./76/LHC.76.RE1111${FileTag}.dat ./LHC.76.REID1111${FileTag}.dat
read
#
summer add ./71/LHC.71.DK1L5${FileTag}.dat ./71/LHC.71.DKRE5${FileTag}.dat ./LHC.71.DKNLO5${FileTag}.dat
#summer add ./71/LHC.71.DK1L5${FileTag}.dat ./71/LHC.71.DKRE5.dat ./LHC.71.DKNLO5${FileTag}.dat
read
#
summer add ./72/LHC.72.DK1L5${FileTag}.dat ./72/LHC.72.DKRE5${FileTag}.dat ./LHC.72.DKNLO5${FileTag}.dat
#summer add ./72/LHC.72.DK1L5${FileTag}.dat ./72/LHC.72.DKRE5.dat ./LHC.72.DKNLO5${FileTag}.dat
read
#
summer add ./LHC.71.DKNLO5${FileTag}.dat ./LHC.72.DKNLO5${FileTag}.dat ./LHC.DKNLO5${FileTag}.dat
read
#
#
#
# alpha=high
summer add ./73/LHC.73.ID2222${FileTag}.dat ./73/LHC.73.RE2222${FileTag}.dat ./LHC.73.REID2222${FileTag}.dat
read
#
summer add ./74/LHC.74.ID2222${FileTag}.dat ./74/LHC.74.RE2222${FileTag}.dat ./LHC.74.REID2222${FileTag}.dat
read
#
summer add ./75/LHC.75.ID2222${FileTag}.dat ./75/LHC.75.RE2222${FileTag}.dat ./LHC.75.REID2222${FileTag}.dat
read
#
summer add ./76/LHC.76.ID2222${FileTag}.dat ./76/LHC.76.RE2222${FileTag}.dat ./LHC.76.REID2222${FileTag}.dat
read
#
summer add ./71/LHC.71.DK1L1${FileTag}.dat ./71/LHC.71.DKRE1${FileTag}.dat ./LHC.71.DKNLO1${FileTag}.dat
#summer add ./71/LHC.71.DK1L1${FileTag}.dat ./71/LHC.71.DKRE1.dat ./LHC.71.DKNLO1${FileTag}.dat
read
#
summer add ./72/LHC.72.DK1L1${FileTag}.dat ./72/LHC.72.DKRE1${FileTag}.dat ./LHC.72.DKNLO1${FileTag}.dat
#summer add ./72/LHC.72.DK1L1${FileTag}.dat ./72/LHC.72.DKRE1.dat ./LHC.72.DKNLO1${FileTag}.dat
read
#
summer add ./LHC.71.DKNLO1${FileTag}.dat ./LHC.72.DKNLO1${FileTag}.dat ./LHC.DKNLO1${FileTag}.dat
read
#
#
#
#
# alpha=low
summer add  ./LHC.1L${FileTag}.dat ./LHC.73.REID1111${FileTag}.dat ./LHC.74.REID1111${FileTag}.dat ./LHC.75.REID1111${FileTag}.dat ./LHC.76.REID1111${FileTag}.dat  ./LHC.NLO1${FileTag}.dat
read
#
#
# alpha=high
summer add ./LHC.1L${FileTag}.dat ./LHC.73.REID2222${FileTag}.dat ./LHC.74.REID2222${FileTag}.dat ./LHC.75.REID2222${FileTag}.dat ./LHC.76.REID2222${FileTag}.dat  ./LHC.NLO2${FileTag}.dat
read
#
#
#
#
#
#
#
#
# putting everything together
echo putting everything into LHC.NLO1all${FileTag}.dat
echo "#LOLO NLO1" > LHC.NLO1all${FileTag}.dat
myjoin ./LHC.LOLO${FileTag}.dat ./LHC.NLO1${FileTag}.dat  >> LHC.NLO1all${FileTag}.dat
read
#
#
echo putting everything into LHC.NLO2all${FileTag}.dat
echo "#LOLO NLO2" > LHC.NLO2all${FileTag}.dat
myjoin ./LHC.LOLO${FileTag}.dat ./LHC.NLO2${FileTag}.dat  >> LHC.NLO2all${FileTag}.dat
read
#
#
#
#
#
#
