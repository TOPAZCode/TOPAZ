// gcc -I /home/schulze/lib/MT2Lib/include/oxbridgekinetics-1.0/  -c calcMT2.cpp

#include "Mt2/ChengHanBisect_Mt2_332_Calculator.h"
// #include <stdlib.h>


extern "C"  {
    void  calcmt2_(double*, double*, double*, double*);
}




void  calcmt2_(double *pA, double *pB, double *pMiss, double *result){
  
    Mt2::ChengHanBisect_Mt2_332_Calculator mt2Calculator;
//     double massOfSystemA =  100;
//     double pxOfSystemA   =  410;
//     double pyOfSystemA   =   20;
//     double massOfSystemB =  150;
//     double pxOfSystemB   = -210;
//     double pyOfSystemB   = -300;
//     double pxMiss        = -200;
//     double pyMiss        =  280;
//     double invis_mass    = 100;
//     Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pxOfSystemA, pyOfSystemA),massOfSystemA);
//     Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pxOfSystemB, pyOfSystemB),massOfSystemB);
//     This should return mT2 = 412.62883811
    
    
    Mt2::LorentzTransverseVector  vis_A(Mt2::TwoVector(pA[0], pA[1]),pA[2]);
    Mt2::LorentzTransverseVector  vis_B(Mt2::TwoVector(pB[0], pB[1]),pB[2]);
    Mt2::TwoVector                pT_Miss(pMiss[0],pMiss[1]);
    
//     mt2Calculator.setDebug(true);
    
    const double mt2 = mt2Calculator.mt2_332(vis_A, vis_B, pT_Miss, pMiss[2]);
    *result = mt2;
    
    
/*    std::cout << "ANSWER: mt2 = " << mt2 
              << " for " << mt2Calculator.algorithmName() << " algorithm"
              << std::endl;  */   
    
    
};


