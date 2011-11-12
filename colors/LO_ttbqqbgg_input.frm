
Global Proj = (SUND(Col5,ColInt15)*SUND(ColInt5,Col15) - 1/NCol*SUND(Col5,Col15)*SUND(ColInt5,ColInt15)) * (SUND(Col6,ColInt16)*SUND(ColInt6,Col16) - 1/NCol*SUND(Col6,Col16)*SUND(ColInt6,ColInt16));

Global [1,Col] = SUND(Col2,Col13)*SUND(Col4,ColInt15)*SUND(ColInt5,ColInt16)*SUND(ColInt6,Col11)   *    Proj;

Global [2,Col] = SUND(Col2,Col13)*SUND(Col4,ColInt16)*SUND(ColInt6,ColInt15)*SUND(ColInt5,Col11)   *    Proj;

Global [3,Col] = SUND(Col2,ColInt15)*SUND(ColInt5,ColInt16)*SUND(ColInt6,Col13)*SUND(Col4,Col11)   *    Proj;

Global [4,Col] = SUND(Col2,ColInt16)*SUND(ColInt6,ColInt15)*SUND(ColInt5,Col13)*SUND(Col4,Col11)   *    Proj;

Global [5,Col] = SUND(Col2,ColInt15)*SUND(ColInt5,Col13)*SUND(Col4,ColInt16)*SUND(ColInt6,Col11)   *    Proj;

Global [6,Col] = SUND(Col2,ColInt16)*SUND(ColInt6,Col13)*SUND(Col4,ColInt15)*SUND(ColInt5,Col11)   *    Proj;

Global [7,Col] = -1/NCol*SUND(Col2,ColInt16)*SUND(ColInt6,ColInt15)*SUND(ColInt5,Col11)*SUND(Col4,Col13)   *    Proj;

Global [8,Col] = -1/NCol*SUND(Col2,ColInt15)*SUND(ColInt5,ColInt16)*SUND(ColInt6,Col11)*SUND(Col4,Col13)   *    Proj;

Global [9,Col] = -1/NCol*SUND(Col2,ColInt15)*SUND(ColInt5,Col11)*SUND(Col4,ColInt16)*SUND(ColInt6,Col13)   *    Proj;

Global [10,Col]= -1/NCol*SUND(Col2,ColInt16)*SUND(ColInt6,Col11)*SUND(Col4,ColInt15)*SUND(ColInt5,Col13)   *    Proj;

Global [11,Col]= -1/NCol*SUND(Col4,ColInt16)*SUND(ColInt6,ColInt15)*SUND(ColInt5,Col13)*SUND(Col2,Col11)   *    Proj;

Global [12,Col]= -1/NCol*SUND(Col4,ColInt15)*SUND(ColInt5,ColInt16)*SUND(ColInt6,Col13)*SUND(Col2,Col11)   *    Proj;





















