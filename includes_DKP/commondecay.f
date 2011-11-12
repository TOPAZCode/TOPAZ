      double complex mom_dec(3,4), hel_dec(3,4), hel_t(1:4), mom_t(1:4)
      double complex zero, ne	
      double complex test_ctm(3,4), test_ctz(3,4), 
     &test_v(3,4), test_r(3,4)
      double precision mt,en,mb,ms,mw, Cf, nc,setalpha, GS
      double precision pii,alphas,el, Gw
      double precision Qqu,qqd
C      parameter(Gw =1.474419561548971d0)
C      parameter(el =0.3028619040941380d0)             
C      parameter(alphaS=0.13d0)
      parameter(pii=3.141592653589793d0)
      parameter(ne=(0d0,1d0))
      parameter(zero=(0d0,0d0))	
C      parameter(qqu= 0.6666666666666666d0)
C      parameter(qqd=-0.3333333333333333d0)
C      parameter(CF = 1.3333333333333333d0)
C      parameter (nc =3.d0)
      integer eve, i_eve, c_eve
      common/kinem/ mom_dec, hel_dec, hel_t, mom_t
      COMMON /mass_color/mt, mb, Mw,setalpha
      common/count/ eve, i_eve,c_eve
      common/testus/test_ctm, test_ctz, test_v, test_r
      common/parameters/Qqd, Qqu, alphas, el, Gw, CF, nc, GS
