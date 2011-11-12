      subroutine numTopDecVirt(Diags)
      implicit none
      include 'VarTopVirt.f'
      include 'commondecay.f'
      include 'intcommon.f'
      double precision wrat, dummy, test
      double complex p1t(4),k2w(4),k1b(4), k3p(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)

      double complex k1bp1t, k2wp1t, k1bk2w
      double complex EpsWk1b, EpsWp1t, EpsWEpsP
      double complex EpsPk1b, EpsPp1t, EpsPk2w
      double complex Diags(8),scf_DKP1L
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
      double complex s21,s22,s23,s24,s25,s26,s27,s28,s29,s30
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1t(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1b(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1b = scf_DKP1L(Dv,EpsW,k1b)
      EpsWp1t = scf_DKP1L(Dv,EpsW,p1t)
      EpsWEpsP = scf_DKP1L(Dv,EpsW,EpsP)
      EpsPk1b = scf_DKP1L(Dv,EpsP,k1b)
      EpsPp1t = scf_DKP1L(Dv,EpsP,p1t)
      EpsPk2w = scf_DKP1L(Dv,EpsP,k2w)
      k1bp1t = scf_DKP1L(Dv,k1b,p1t)
      k1bk2w = scf_DKP1L(Dv,k1b,k2w)
      k2wp1t = scf_DKP1L(Dv,k2w,p1t)
!      write(*,*) "Top-Decay-Virtual-Photon switched off here"
      include 'TopVirtuals.f'
!      stop
      end subroutine numTopDecVirt


      subroutine numATopDecVirt(Diags)
      implicit none
      include 'VarATopVirt.f'
      include 'commondecay.f'
      include 'intcommonA.f'
      double complex p1tb(4),k2w(4),k1bb(4)
      double complex EpsP(4),EpsW(4),SpB(4),SpT(4)

      double complex k1bbp1tb, k2wp1tb, k1bbk2w
      double complex EpsWk1bb, EpsWp1tb, EpsWEpsP
      double complex EpsPk1bb, EpsPp1tb, EpsPk2w
      double complex Diags(8),scf_DKP1L
      double complex s1,s2,s3,s4,s5,s6,s7,s8,s9,s10
      double complex s11,s12,s13,s14,s15,s16,s17,s18,s19,s20
      double complex s21,s22,s23,s24,s25,s26,s27,s28,s29,s30
      double complex s31,s32,s33,s34,s35,s36,s37,s38,s39,s40
      integer i, Dv
      Dv = 4
      do i = 1,4
         p1tb(i)=mom_t(i)
         k2w(i)=mom_dec(2,i)
         k1bb(i)=mom_dec(1,i)
         EpsP(i)=hel_dec(3,i)
         EpsW(i)=hel_dec(2,i)
         SpB(i)=hel_dec(1,i)
      enddo
      EpsWk1bb = scf_DKP1L(Dv,EpsW,k1bb)
      EpsWp1tb = scf_DKP1L(Dv,EpsW,p1tb)
      EpsWEpsP = scf_DKP1L(Dv,EpsW,EpsP)
      EpsPk1bb = scf_DKP1L(Dv,EpsP,k1bb)
      EpsPp1tb = scf_DKP1L(Dv,EpsP,p1tb)
      EpsPk2w = scf_DKP1L(Dv,EpsP,k2w)
      k1bbp1tb = scf_DKP1L(Dv,k1bb,p1tb)
      k1bbk2w = scf_DKP1L(Dv,k1bb,k2w)
      k2wp1tb = scf_DKP1L(Dv,k2w,p1tb)
!      write(*,*) "AntiTop-Decay-Virtual-Photon switched off here"
!      stop
      include 'ATopVirtuals.f'
      end subroutine numATopDecVirt


      double complex function scf_DKP1L(Dv,ap1,ap2)
        implicit none
        integer Dv
        double complex ap1(Dv),ap2(Dv)
        double complex r1

        call sc_DKP1L(Dv,ap1,ap2,r1)
        scf_DKP1L = r1
        return
      end function

      subroutine sc_DKP1L(n,x,y,r)
        implicit none
        integer i,n
        double complex x(*),y(*)
        double complex r

        r=dcmplx(0d0,0d0)

        do i=1, n
           if (i.eq.1) then
              r = r + x(i)*y(i)
           else
              r = r - x(i)*y(i)
           endif
        enddo

        return
      end subroutine


