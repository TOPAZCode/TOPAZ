*  generate color correlation matrix

#define WorkPath "/home/schulze/projects/ttbjets/colors/"
#define LibPath "/home/schulze/lib/FeynArtsToForm/"
*_____________________________________________________________________________________________*
#define InputColorFile "/home/schulze/projects/ttbjets/colors/LO_ttbgggg_input.frm"
#define OutputFile "/home/schulze/projects/ttbjets/colors/LO_ttbgggg_output"
*_____________________________________________________________________________________________*
#define setDSTm4 "0"
#define setDST   "4"

#define InterfereDiags1 "{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}";
#define InterfereDiags2 "{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24}";
*_____________________________________________________________________________________________*
#include `LibPath'header2.frm;
#include `InputColorFile';
#include `LibPath'QCDColorInterf.frm;
*_____________________________________________________________________________________________*

id NCol=3;
id NCol**(-2)=1/9;
id NCol**(-1)=1/3;


Format doublefortran;
Print;
.sort

#do TheCcAmp = `InterfereDiags2'
#do TheAmp = `InterfereDiags1'
  #write <`OutputFile'.dat> "ColLO_ttbgggg[`TheAmp',{`TheCcAmp'}] = (%E); ", [`TheAmp',`TheCcAmp',Col]
#enddo
#enddo
#write <`OutputFile'.dat> "\n"

.store
Off statistics;

#system sed "s/ 0/ 0.d0/g" `OutputFile'.dat | sed "s/\[/\(/g" | sed "s/\]/\)/g" | sed "s/;/ /g" > `OutputFile'.f90;
#system rm `OutputFile'.dat;

.end
