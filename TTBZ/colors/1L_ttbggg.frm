*  generate color correlation matrix

#define WorkPath "/home/schulze/projects/ttbjets/colors/"
#define LibPath "/home/schulze/lib/FeynArtsToForm/"
*_____________________________________________________________________________________________*
#define InputColorFile "/home/schulze/projects/ttbjets/colors/1L_ttbggg_input.frm"
#define OutputFile "/home/schulze/projects/ttbjets/colors/1L_ttbggg_output"
*_____________________________________________________________________________________________*
#define setDSTm4 "0"
#define setDST   "4"

#define InterfereDiags1 "{1,2,3,4,5,6}";
#define InterfereDiags2 "{10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44}";
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

#do TheAmp = `InterfereDiags1'
#do TheCcAmp = `InterfereDiags2'
  #write <`OutputFile'.dat> "Col1L_ttbggg[`TheAmp',{`TheCcAmp'-9}] = (%E); ", [`TheAmp',`TheCcAmp',Col]
#enddo
#enddo
#write <`OutputFile'.dat> "\n"

.store
Off statistics;

#system sed "s/ 0/ 0.d0/g" `OutputFile'.dat | sed "s/\[/\(/g" | sed "s/\]/\)/g" | sed "s/;/ /g" > `OutputFile'.f90;
#system rm `OutputFile'.dat;

.end
