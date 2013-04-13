*  generate color correlation matrix

#define WorkPath "/home/schulze/projects/TOPAZ/colors/"
#define LibPath "/home/schulze/lib/FeynArtsToForm/"
*_____________________________________________________________________________________________*
#define InputColorFile "/home/schulze/projects/TOPAZ/colors/Dipole_tDKgg_input.frm"
#define OutputFile "/home/schulze/projects/TOPAZ/colors/Dipole_tDKgg_output"
*_____________________________________________________________________________________________*
#define setDSTm4 "0"
#define setDST   "4"

#define InterfereDiags1 "1,1";
#define InterfereDiags2 "1,1";


*----  T_t * T_b
#define ColorCorrelator "SUNT(GluInt1,ColP2,Col2)*(-SUNT(GluInt1,Col1,ColP1))*SUND(Glu1,GluP1)";

*----  T_t * T_g
*#define ColorCorrelator "-SUNT(Glu9,Col1,ColP1)*i_*SUNF(GluP1,Glu9,Glu1)*SUND(Col2,ColP2)";

*----  T_b * T_g
*#define ColorCorrelator "SUNT(GluInt1,ColP2,Col2)*i_*SUNF(GluP1,GluInt1,Glu1)*SUND(Col1,ColP1)";



*____________________________________________________________________________________________*
#include `LibPath'header2.frm;
#include `InputColorFile';
#include `LibPath'QCDColorInterf.frm;
*_____________________________________________________________________________________________*

id NCol=3;
id NCol**(-2)=1/9;
id NCol**(-1)=1/3;


* unknown correction factor
multiply 1/4;




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
