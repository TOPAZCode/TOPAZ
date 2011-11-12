Global [1,Col] = SUNT(Glu3,Col2,ColInt1) * SUNT(Glu4, ColInt1, Col1);
Global [2,Col] = SUNT(Glu4,Col2,ColInt1) * SUNT(Glu3, ColInt1, Col1);

Global [10,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt3))*(i_*SUNF(GluInt3,Glu3,GluInt2));
Global [11,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt3))*(i_*SUNF(GluInt3,Glu4,GluInt2));

Global [12,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt2));
Global [13,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt2));

Global [14,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(Glu3,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * SUND(GluInt1,GluInt2);
Global [15,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(Glu4,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * SUND(GluInt1,GluInt2);


* fermion loop contributions (to be multiplied by (-1)*N_f):
Global [16,Col] = SUNT(Glu3,Col2,ColInt1) * SUNT(Glu4,ColInt1,Col1);
Global [17,Col] = SUNT(Glu4,Col2,ColInt1) * SUNT(Glu3,ColInt1,Col1);

Global [18,Col] = 1/NCol * SUNT(Glu3,ColInt2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2) * SUND(Col2,Col1);


