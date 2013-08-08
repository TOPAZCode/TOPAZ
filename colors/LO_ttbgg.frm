*  generate color correlation matrix

#define WorkPath "/home/schulze/projects/TOPAZ/colors"
#define LibPath "/home/schulze/lib/FeynArtsToForm/"
*_____________________________________________________________________________________________*
#define InputColorFile "/home/schulze/projects/TOPAZ/colors/LO_ttbgg_input.frm"
#define OutputFile "/home/schulze/projects/TOPAZ/colors/LO_ttbgg_output"
*_____________________________________________________________________________________________*
#define setDSTm4 "0"
#define setDST   "4"

#define InterfereDiags1 "{1,2}";
#define InterfereDiags2 "{1,2}";



*#define ColorCorrelator "i_*SUND(Col2, ColP2)*SUND(Glu4, GluP4)*SUNF(GluInt20, Glu3, GluP3)*SUNT(GluInt20, ColP1, Col1)"
*#define ColorCorrelator "i_*SUND(Col2, ColP2)*SUND(Glu4, GluP4)*SUNF(GluInt20, Glu3, GluP3)*SUNT(GluInt20, ColP1, Col1)"


*____________________________________________________________________________________________*
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
  #write <`OutputFile'.dat> "ColLO_ttbgg[`TheAmp',{`TheCcAmp'}] = (%E); ", [`TheAmp',`TheCcAmp',Col]
#enddo
#enddo
#write <`OutputFile'.dat> "\n"

.store
Off statistics;

#system sed "s/ 0/ 0.d0/g" `OutputFile'.dat | sed "s/\[/\(/g" | sed "s/\]/\)/g" | sed "s/;/ /g" > `OutputFile'.f90;
#system rm `OutputFile'.dat;

.end
