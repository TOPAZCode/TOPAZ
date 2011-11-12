      coeff5(CutNum,1)= PentCuts%Coeff(CutNum,0)
      mass5(CutNum,1) = IntPart( PentCuts%CutProp(CutNum,1) )%Mass
      mass5(CutNum,2) = IntPart( PentCuts%CutProp(CutNum,2) )%Mass
      mass5(CutNum,3) = IntPart( PentCuts%CutProp(CutNum,3) )%Mass
      mass5(CutNum,4) = IntPart( PentCuts%CutProp(CutNum,4) )%Mass
      mass5(CutNum,5) = IntPart( PentCuts%CutProp(CutNum,5) )%Mass
      propv5(CutNum,1:4)  = (0d0,0d0)
      propv5(CutNum,5:8)  = KMom(1,1:4)
      propv5(CutNum,9:12) = KMom(2,1:4)
      propv5(CutNum,13:16)= KMom(3,1:4)
      propv5(CutNum,17:20)= KMom(4,1:4)
!-----------------------------------------------------------------------
      coeff4(CutNum,1:5) = QuadCuts%Coeff(CutNum,0:4)
      refvect4(CutNum,1:4) = NMom(1,1:4)
      mass4(CutNum,1) = IntPart( QuadCuts%CutProp(CutNum,1) )%Mass
      mass4(CutNum,2) = IntPart( QuadCuts%CutProp(CutNum,2) )%Mass
      mass4(CutNum,3) = IntPart( QuadCuts%CutProp(CutNum,3) )%Mass
      mass4(CutNum,4) = IntPart( QuadCuts%CutProp(CutNum,4) )%Mass
      propv4(CutNum,1:4)  = (0d0,0d0)
      propv4(CutNum,5:8)  = KMom(1,1:4)
      propv4(CutNum,9:12) = KMom(2,1:4)
      propv4(CutNum,13:16)= KMom(3,1:4)
!-----------------------------------------------------------------------
      coeff3(CutNum,1) = TripCuts%Coeff(CutNum,0)
      coeff3(CutNum,2) = TripCuts%Coeff(CutNum,1)
      coeff3(CutNum,3) = TripCuts%Coeff(CutNum,2)
      coeff3(CutNum,4) = TripCuts%Coeff(CutNum,3)
      coeff3(CutNum,5) = TripCuts%Coeff(CutNum,4)
      coeff3(CutNum,6) = TripCuts%Coeff(CutNum,5)
      coeff3(CutNum,7) = TripCuts%Coeff(CutNum,6)
      coeff3(CutNum,8) = TripCuts%Coeff(CutNum,7)
      coeff3(CutNum,9) = TripCuts%Coeff(CutNum,8)
      coeff3(CutNum,10)= TripCuts%Coeff(CutNum,9)
      refvect3(CutNum,1:4) = NMom(1,1:4)
      refvect3(CutNum,5:8) = NMom(2,1:4)
      mass3(CutNum,1) = IntPart( TripCuts%CutProp(CutNum,1) )%Mass
      mass3(CutNum,2) = IntPart( TripCuts%CutProp(CutNum,2) )%Mass
      mass3(CutNum,3) = IntPart( TripCuts%CutProp(CutNum,3) )%Mass
      propv3(CutNum,1:4)  = (0d0,0d0)
      propv3(CutNum,5:8)  = KMom(1,1:4)
      propv3(CutNum,9:12) = KMom(2,1:4)
!-----------------------------------------------------------------------
         coeff2(CutNum,1) = DoubCuts%Coeff(CutNum,0)
         coeff2(CutNum,2) = DoubCuts%Coeff(CutNum,1)
         coeff2(CutNum,3) = DoubCuts%Coeff(CutNum,2)
         coeff2(CutNum,4) = DoubCuts%Coeff(CutNum,3)
         coeff2(CutNum,5) = DoubCuts%Coeff(CutNum,4)
         coeff2(CutNum,6) = DoubCuts%Coeff(CutNum,5)
         coeff2(CutNum,7) = DoubCuts%Coeff(CutNum,6)
         coeff2(CutNum,8) = DoubCuts%Coeff(CutNum,7)
         coeff2(CutNum,9) = DoubCuts%Coeff(CutNum,8)
         coeff2(CutNum,10)= DoubCuts%Coeff(CutNum,9)
         refvect2(CutNum,1:4)  = NMom(1,1:4)
         refvect2(CutNum,5:8)  = NMom(2,1:4)
         refvect2(CutNum,9:12) = NMom(3,1:4)
         mass2(CutNum,1) = IntPart( DoubCuts%CutProp(CutNum,1) )%Mass
         mass2(CutNum,2) = IntPart( DoubCuts%CutProp(CutNum,2) )%Mass
         propv2(CutNum,1:4)  = (0d0,0d0)
         propv2(CutNum,5:8)  = KMom(1,1:4)
!-----------------------------------------------------------------------
         coeff1(CutNum,1) = SingCuts%Coeff(CutNum,0)
         mass1(CutNum,1)  = IntPart( SingCuts%CutProp(CutNum,1) )%Mass





