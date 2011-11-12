Global [1,Col] = SUNT(Glu5,Col4,Col1) * SUND(Col2,Col3);
Global [2,Col] = 1/NCol * SUNT(Glu5,Col2,Col1) * SUND(Col4,Col3);
Global [3,Col] = SUNT(Glu5,Col2,Col3) * SUND(Col1,Col4);
Global [4,Col] = 1/NCol * SUNT(Glu5,Col4,Col3) * SUND(Col2,Col1);

Global [10,Col] = SUNT(Glu5,Col4,Col1) * SUND(Col2,Col3);
Global [11,Col] = SUNT(Glu5,Col2,Col1) * SUND(Col4,Col3);
Global [12,Col] = SUNT(Glu5,Col2,Col3) * SUND(Col1,Col4);
Global [13,Col] = SUNT(Glu5,Col4,Col3) * SUND(Col2,Col1);

