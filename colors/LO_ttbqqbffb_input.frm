Global [1,Col] = SUND(Col2,Col3)*SUND(Col4,Col5)*SUND(Col6,Col1);
Global [2,Col] = SUND(Col2,Col5)*SUND(Col6,Col3)*SUND(Col4,Col1);
Global [3,Col] = -1/NCol*( SUND(Col4,Col3)*SUND(Col2,Col5)*SUND(Col6,Col1) + SUND(Col2,Col1)*SUND(Col6,Col3)*SUND(Col4,Col5) - 1/NCol*SUND(Col2,Col1)*SUND(Col6,Col5)*SUND(Col4,Col3) );

Global [4,Col] = -1/NCol*( SUND(Col5,Col6)*SUND(Col2,Col3)*SUND(Col4,Col1) + SUND(Col3,Col4)*SUND(Col6,Col1)*SUND(Col2,Col5) - 1/NCol*SUND(Col2,Col1)*SUND(Col6,Col5)*SUND(Col4,Col3) );


Global [5,Col] = -1/NCol*( SUND(Col5,Col6)*SUND(Col2,Col3)*SUND(Col4,Col1) + SUND(Col2,Col1)*SUND(Col6,Col3)*SUND(Col4,Col5) - 1/NCol*SUND(Col2,Col1)*SUND(Col6,Col5)*SUND(Col4,Col3) );

