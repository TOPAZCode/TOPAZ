
! ! ! ! ! ! ! ! ! ! ! !    C  H  E C  K   ! ! ! ! ! ! ! ! ! ! ! !

call coupsm(0)
if( NLOParam.ne.0 ) call Error("NLOParam has to be zero, i.e. alpha_s=0.13!")
if( PDFSet.ne.2 ) call Error("call with CTEQ PDFs for alpha_s(LO)=0.13")

IF( Contr.eq.1 ) THEN


      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do GluHel2=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,GluHel2,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( ATop(GluHel1,GluHel2,2)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) + ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0)) + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) + ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact.
      call STB_BBEMVEBGG((/MomExt(1:4,3)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause






ELSEIF( Contr.eq.2 ) THEN


      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do GluHel2=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,GluHel2,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( ATop(GluHel1,GluHel2,2)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) + ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0)) + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) + ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
!       NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact.! already included in Andreas' code
      call STB_BBUBDGG((/MomExt(1:4,3)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause








ELSEIF( Contr.eq.8 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do GluHel2=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,GluHel2,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( Top(GluHel1,GluHel2,2)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) + ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0)) + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) + ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact.
      call ST_BEPVEGG((/MomExt(1:4,4)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause






ELSEIF( Contr.eq.9 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_
      do GluHel2=plus_,minus_
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,GluHel2,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( Top(GluHel1,GluHel2,2)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) + ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0)) + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) + ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
!       NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact.! already included in Andreas' code
      call ST_BDBUGG((/MomExt(1:4,4)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause










ELSEIF( Contr.eq.11 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1) * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBBUUB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause






ELSEIF( Contr.eq.12 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( ATop(GluHel1,0,2)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) +  ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0))
                    if(GluHel1.eq.plus_) NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) +  ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact
      call STB_EMVEBBBBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,8)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause





ELSEIF( Contr.eq.13 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1) * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEBUUB((/MomExt(1:4,4)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,11)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause




ELSEIF( Contr.eq.14 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( Top(GluHel1,0,2)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) +  ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0))
                    if(GluHel1.eq.plus_) NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) +  ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact
      call ST_EPVEBBBB((/MomExt(1:4,4)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,11)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause






ELSEIF( Contr.eq.15 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1) * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_BBUBDCCB((/MomExt(1:4,3)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause





ELSEIF( Contr.eq.16 ) THEN
      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( ATop(GluHel1,0,2)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) +  ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0))
                    if(GluHel1.eq.plus_) NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) +  ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact
      call STB_BBUBDUUB((/MomExt(1:4,3)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8

      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause







ELSEIF( Contr.eq.17 ) THEN
      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( ATop(GluHel1,0,2)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) +  ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0))
                    if(GluHel1.eq.plus_) NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) +  ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm.fact
      call STB_BBUBDDDB((/MomExt(1:4,3)*100d0,MomExt(1:4,5)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8

      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_EPVEB((/MomExt(1:4,4)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0,MomExt(1:4,10)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause






ELSEIF( Contr.eq.18 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""




      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1) * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call ST_BDBUCCB((/MomExt(1:4,4)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause





ELSEIF( Contr.eq.19 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""



      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( Top(GluHel1,0,2)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) +  ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0))
                    if(GluHel1.eq.plus_) NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) +  ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm fact.
      call ST_BDBUDDB((/MomExt(1:4,4)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause










ELSEIF( Contr.eq.20 ) THEN

      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call vbarSpi(ExtParticle(1)%Mom(1:4),m_Top,iHel,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = psp1_( ATop(GluHel1,0,1)%Pol(1:4) ,ExtParticle(1)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + 1d0 * TreeResult(1,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      call STB_EMVEBBB((/MomExt(1:4,3)*100d0,MomExt(1:4,6)*100d0,MomExt(1:4,7)*100d0,MomExt(1:4,5)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d0
      print *, "ATop Spinor:",NLO_Res_Unpol
      print *, "MadGraph:   ",dummy(5)
      print *, "ratio:      ",NLO_Res_Unpol/dummy(5)
      print *, ""



      NLO_Res_Unpol = (0d0,0d0)
      do GluHel1=plus_,minus_      ! loop over additional q-qbar polarization
      do iHel=1,-1,-2! top pol
                    call uSpi(ExtParticle(2)%Mom(1:4),m_Top,iHel,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = psp1_( Top(GluHel1,0,1)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(1,0) = TreeResult(1,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(1,0) = TreeResult(1,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    TreeResult(2,0) = psp1_( Top(GluHel1,0,2)%Pol(1:4) ,ExtParticle(2)%Pol(1:4))
                    TreeResult(2,0) = TreeResult(2,0)/(2d0*m_Top)! remove this factor coming from (pslash-m_top) projection
                    TreeResult(2,0) = TreeResult(2,0) *dsqrt(2d0*Ga_Top(0)*m_Top) *dsqrt(2d0*Ga_W*m_W) /(Ga_W*m_W)!  remove NWAFactor_Top, remove WProp(NWA), add full prop. for pW2=MW2

                    NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,1)*TreeResult(1,0)*dconjg(TreeResult(1,0)) +  ColCorrDK(2,2)*TreeResult(2,0)*dconjg(TreeResult(2,0))
                    if(GluHel1.eq.plus_) NLO_Res_UnPol = NLO_Res_UnPol + ColCorrDK(1,2)*TreeResult(1,0)*dconjg(TreeResult(2,0)) +  ColCorrDK(2,1)*TreeResult(2,0)*dconjg(TreeResult(1,0))
      enddo!helicity loop
      enddo!helicity loop
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! spin avg.
      NLO_Res_Unpol = NLO_Res_Unpol  * 1d0/2d0! symm fact.
      call ST_BDBUUUB((/MomExt(1:4,4)*100d0,MomExt(1:4,8)*100d0,MomExt(1:4,9)*100d0,MomExt(1:4,10)*100d0,MomExt(1:4,11)*100d0,MomExt(1:4,12)*100d0/),dummy(5))
      dummy(5) = dummy(5)*1d8
      print *, "Top Spinor:",NLO_Res_Unpol
      print *, "MadGraph:  ",dummy(5)
      print *, "ratio:     ",NLO_Res_Unpol/dummy(5)
      print *, ""

      pause














ENDIF