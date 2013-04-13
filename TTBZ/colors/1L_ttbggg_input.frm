Global [1,Col] = SUNT(Glu3,Col2,ColInt1) * SUNT(Glu4,ColInt1,ColInt2) * SUNT(Glu5,ColInt2,Col1);
Global [2,Col] = SUNT(Glu3,Col2,ColInt1) * SUNT(Glu5,ColInt1,ColInt2) * SUNT(Glu4,ColInt2,Col1);
Global [3,Col] = SUNT(Glu4,Col2,ColInt1) * SUNT(Glu3,ColInt1,ColInt2) * SUNT(Glu5,ColInt2,Col1);
Global [4,Col] = SUNT(Glu4,Col2,ColInt1) * SUNT(Glu5,ColInt1,ColInt2) * SUNT(Glu3,ColInt2,Col1);
Global [5,Col] = SUNT(Glu5,Col2,ColInt1) * SUNT(Glu3,ColInt1,ColInt2) * SUNT(Glu4,ColInt2,Col1);
Global [6,Col] = SUNT(Glu5,Col2,ColInt1) * SUNT(Glu4,ColInt1,ColInt2) * SUNT(Glu3,ColInt2,Col1);

Global [10,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu5,GluInt3))*(i_*SUNF(GluInt3,Glu4,GluInt4))*(i_*SUNF(GluInt4,Glu3,GluInt2));
Global [11,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt3))*(i_*SUNF(GluInt3,Glu5,GluInt4))*(i_*SUNF(GluInt4,Glu3,GluInt2));
Global [12,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu5,GluInt3))*(i_*SUNF(GluInt3,Glu3,GluInt4))*(i_*SUNF(GluInt4,Glu4,GluInt2));
Global [13,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt3))*(i_*SUNF(GluInt3,Glu5,GluInt4))*(i_*SUNF(GluInt4,Glu4,GluInt2));
Global [14,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt3))*(i_*SUNF(GluInt3,Glu3,GluInt4))*(i_*SUNF(GluInt4,Glu5,GluInt2));
Global [15,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(GluInt1,ColInt1,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt3))*(i_*SUNF(GluInt3,Glu4,GluInt4))*(i_*SUNF(GluInt4,Glu5,GluInt2));

Global [16,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu5,GluInt3))*(i_*SUNF(GluInt3,Glu4,GluInt2));
Global [17,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt3))*(i_*SUNF(GluInt3,Glu5,GluInt2));
Global [18,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu5,GluInt3))*(i_*SUNF(GluInt3,Glu3,GluInt2));
Global [19,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt3))*(i_*SUNF(GluInt3,Glu5,GluInt2));
Global [20,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt3))*(i_*SUNF(GluInt3,Glu3,GluInt2));
Global [21,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(GluInt1,ColInt2,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt3))*(i_*SUNF(GluInt3,Glu4,GluInt2));

Global [22,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(Glu3,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * (i_*SUNF(GluInt1,Glu5,GluInt2));
Global [23,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(Glu3,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt2));
Global [24,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(Glu4,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * (i_*SUNF(GluInt1,Glu5,GluInt2));
Global [25,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(Glu4,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt2));
Global [26,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(Glu5,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * (i_*SUNF(GluInt1,Glu4,GluInt2));
Global [27,Col] = SUNT(GluInt2,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(Glu5,ColInt2,ColInt3)*SUNT(GluInt1,ColInt3,Col1) * (i_*SUNF(GluInt1,Glu3,GluInt2));

Global [28,Col] = SUNT(GluInt1,Col2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(Glu4,ColInt2,ColInt3)*SUNT(Glu3,ColInt3,ColInt4)*SUNT(GluInt1,ColInt4,Col1);
Global [29,Col] = SUNT(GluInt1,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(Glu5,ColInt2,ColInt3)*SUNT(Glu3,ColInt3,ColInt4)*SUNT(GluInt1,ColInt4,Col1);
Global [30,Col] = SUNT(GluInt1,Col2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(Glu3,ColInt2,ColInt3)*SUNT(Glu4,ColInt3,ColInt4)*SUNT(GluInt1,ColInt4,Col1);
Global [31,Col] = SUNT(GluInt1,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(Glu5,ColInt2,ColInt3)*SUNT(Glu4,ColInt3,ColInt4)*SUNT(GluInt1,ColInt4,Col1);
Global [32,Col] = SUNT(GluInt1,Col2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(Glu3,ColInt2,ColInt3)*SUNT(Glu5,ColInt3,ColInt4)*SUNT(GluInt1,ColInt4,Col1);
Global [33,Col] = SUNT(GluInt1,Col2,ColInt1)*SUNT(Glu3,ColInt1,ColInt2)*SUNT(Glu4,ColInt2,ColInt3)*SUNT(Glu5,ColInt3,ColInt4)*SUNT(GluInt1,ColInt4,Col1);


* fermion loop contributions (to be multiplied by (-1)*N_f):
Global [34,Col] = SUNT(Glu3,Col2,ColInt1) * SUNT(Glu4,ColInt1,ColInt2) * SUNT(Glu5,ColInt2,Col1);
Global [35,Col] = SUNT(Glu3,Col2,ColInt1) * SUNT(Glu5,ColInt1,ColInt2) * SUNT(Glu4,ColInt2,Col1);
Global [36,Col] = SUNT(Glu4,Col2,ColInt1) * SUNT(Glu3,ColInt1,ColInt2) * SUNT(Glu5,ColInt2,Col1);
Global [37,Col] = SUNT(Glu4,Col2,ColInt1) * SUNT(Glu5,ColInt1,ColInt2) * SUNT(Glu3,ColInt2,Col1);
Global [38,Col] = SUNT(Glu5,Col2,ColInt1) * SUNT(Glu3,ColInt1,ColInt2) * SUNT(Glu4,ColInt2,Col1);
Global [39,Col] = SUNT(Glu5,Col2,ColInt1) * SUNT(Glu4,ColInt1,ColInt2) * SUNT(Glu3,ColInt2,Col1);



Global [40,Col] = 1/NCol * SUNT(Glu3,ColInt2,ColInt1)*SUNT(Glu4,ColInt1,ColInt2) * SUNT(Glu5,Col2,Col1);
Global [41,Col] = 1/NCol * SUNT(Glu3,ColInt2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2) * SUNT(Glu4,Col2,Col1);
Global [42,Col] = 1/NCol * SUNT(Glu4,ColInt2,ColInt1)*SUNT(Glu5,ColInt1,ColInt2) * SUNT(Glu3,Col2,Col1);

Global [43,Col] = 1/NCol * SUNT(Glu3,ColInt3,ColInt1)*SUNT(Glu4,ColInt1,ColInt2)*SUNT(Glu5,ColInt2,ColInt3) * SUND(Col2,Col1);
Global [44,Col] = 1/NCol * SUNT(Glu3,ColInt3,ColInt1)*SUNT(Glu5,ColInt1,ColInt2)*SUNT(Glu4,ColInt2,ColInt3) * SUND(Col2,Col1);











