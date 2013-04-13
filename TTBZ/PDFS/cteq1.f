      FUNCTION CteqPd (Iset, Iparton, X, Q, Irt)

C     CTEQ distribution function in a parametrized form.  

C     (No data tables are needed.)

C The returned function value is the PROBABILITY density for a given FLAVOR.

C  !! A companion function (next module), which this one depends on, 
C  !!        CtqPdf (Iset, Iparton, X, Q, Irt)
C  !! gives the VALENCE and SEA MOMENTUM FRACTION distributions. 

C  \\  A parallel (independent) program CtqPds (not included in this file) 
C  ||  in Subroutine form is also available. 
C  ||  It returns ALL the parton flavors at once in an array form.
C  //  See details in that separate file if you are interested.

C Ref.: "CTEQ Parton Distributions and Flavor Dependence of the Sea Quarks"
C       by: J. Botts, J.G. Morfin, J.F. Owens, J. Qiu, W.K. Tung & H. Weerts
C       MSUHEP-92-27, Fermilab-Pub-92/371, FSU-HEP-92-1225, ISU-NP-92-17
C       Now published in Phys. Lett.B304, 159 (1993).

C     Since this is an initial distribution, and there may be updates, it is 
C     useful for the authors to maintain a record of the distribution list.
C     Please do not freely distribute this program package; instead, refer any 
C     interested colleague to direct their request for a copy to:
C     Botts@msupa.pa.msu.edu  or  Botts@msupa (bitnet)  or  MSUHEP::Botts

C    If you have any questions concerning these distributions, direct inquires 
C    to Jim Botts or Wu-Ki Tung (username Tung at same E-mail nodes as above).

C$Header: /users/wkt/1hep/0cteq/RCS/CteqPd.f,v 1.2 93/02/26 10:42:43 wkt Exp $
C$Log:	CteqPd.f,v $
c Revision 1.2  93/02/26  10:42:43  wkt
c Version with heavy quark threshold factor and faster algorithm.
c 
c Revision 1.1  93/02/14  17:30:21  botts
c The new Faster version.
c Revision 1.0  93/02/08  18:35:25  wkt
c Initial revision

C   This function returns the CTEQ parton distributions f^Iset_Iprtn/proton
C     where Iset (= 1, 2, ..., 5) is the set label; 

C       Name convention for CTEQ distributions:  CTEQnSx  where
C           n : version number                      (currently n = 1)
C           S : factorization scheme label: = [M D L] for [MS-bar DIS LO]  
c               resp.
C           x : special characteristics, if any 
C                    (e.g. S for singular gluon, L for "LEP lambda value")

C   Iprtn  is the parton label (6, 5, 4, 3, 2, 1, 0, -1, ......, -6)
C                          for (t, b, c, s, d, u, g, u_bar, ..., t_bar)

C X, Q are the usual x, Q; Irt is a return error code (not implemented yet).

C --> Iset = 1, 2, 3, 4, 5 correspond to the following CTEQ global fits:
C     cteq1M, cteq1MS, cteq1ML, cteq1D, cteq1L  respectively.

C --> QCD parameters for parton distribution set Iset can be obtained inside
C         the user's program by:
C     Dum = Prctq1 
C    >        (Iset, Iord, Ischeme, MxFlv,
C    >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
C    >         Xmin, Qini, Qmax, ExpNor)
C     where all but the first argument are output parameters.
C     They should be self-explanary -- see details in next module.

C     The range of (x, Q) used in this round of global analysis is, approxi-
C     mately,  0.01 < x < 0.75 ; and 4 GeV^2 < Q^2 < 400 GeV^2.

C    The range of (x, Q) used in the reparametrization of the QCD evolved
C    parton distributions is 10E-5 < x < 1 ; 2 GeV < Q < 1 TeV.  The  
C    functional form of this parametrization is:

C      A0 * x^A1 * (1-x)^A2 * (1 + A3 * x^A4) * [log(1+1/x)]^A5

C   with the A'coefficients being smooth functions of Q.  For heavy quarks,
C   an additional threshold factor is applied which simulate the Q-dependence
C   of the QCD evolution in that region.

C   Since this function is positive definite and smooth, it provides sensible
C   extrapolations of the parton distributions if they are called beyond
C   the original range in an application. There is no artificial boundaries
C   or sharp cutoff's.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Ifl = Iparton
      JFL = ABS(Ifl)
C                                                             Valence
      IF (Ifl .LE. 0) THEN
        VL = 0
      ELSEIF (Ifl .LE. 2) THEN
        VL = CtqPdf(Iset, Ifl, X, Q, Irt)
      ELSE
        VL = 0
      ENDIF
C                                                                         Sea
      SEA = CtqPdf (Iset, -JFL, X, Q, Irt)
C                                              Full (probability) Distribution 
      CteqPd = (VL + SEA) / X
       
      Return
C                         *************************
      END

      FUNCTION CtqPdf (Iset, Iprtn, XX, QQ, Irt)

C            Returns xf(x,Q) -- the momentum fraction distribution !!
C            Returns valence and sea rather than combined flavor distr.

C            Iset : PDF set label

C            Iprtn  : Parton label:   2, 1 = d_ and u_ valence
C                                     0 = gluon
C                            -1, ... -6 = u, d, s, c, b, t sea quarks

C            XX  : Bjorken-x
C            QQ  : scale parameter "Q"
C            Irt : Return code

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (Nex = 5, MxFl = 6, Npn = 3, Nst = 30, Nexpt=20)
      Parameter (Nst4 = Nst*4)

      DIMENSION
     >   Iord(Nst), Isch(Nst), Nqrk(Nst),Alm(Nst)
     > , Vlm(4:6,Nst), Qms(4:6, Nst)
     > , Xmn(Nst), Qmn(Nst), Qmx(Nst), Nexp(Nexpt)
     > , Mex(Nst), Mpn(Nst), ExpN(Nexpt, Nst), ExpNor(Nexpt)

c                                                             CTEQ1M
      DATA 
     >  Isch(1), Iord(1), Nqrk(1), Alm(1) /  1,  2,  6,  .152 / 
     >  (Vlm(I,1), I=4,6) / .231,    .152,    .059  /
     >  (Qms(I,1), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(1), Qmn(1), Qmx(1) /  1.E-5,  2.00,  1.E3  /
     >  Mex(1), Mpn(1), Nexp(1) /  5, 3, 8  /
     >  (ExpN(I, 1), I=1,8)
     >  / 0.989, 1.00, 1.02, 0.978, 1.10, 0.972, 0.987, 0.846 /
c                                                             CTEQ1MS
      DATA 
     >  Isch(2), Iord(2), Nqrk(2), Alm(2) /  1,  2,  6, .152  / 
     >  (Vlm(I,2), I=4,6) / .231,    .152,    .059  /
     >  (Qms(I,2), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(2), Qmn(2), Qmx(2) /  1.E-5,  2.00,  1.E3  /
     >  Mex(2), Mpn(2), Nexp(2) /  5, 3, 8  /
     >  (ExpN(I, 2), I=1,8 )
     >  / 0.989, 1.00, 1.02, 0.984, 1.05, 0.891, 0.923, 0.824 /
c                                                             CTEQ1ML
      DATA 
     >  Isch(3), Iord(3), Nqrk(3), Alm(3) /  1,  2,  6, .220  / 
     >  (Vlm(I,3), I=4,6) / .322,    .220,     .088  /
     >  (Qms(I,3), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(3), Qmn(3), Qmx(3) /  1.E-5,  2.00,  1.E3  /
     >  Mex(3), Mpn(3), Nexp(3) /  5, 3, 8 /
     >  (ExpN(I, 3), I=1,8 )
     >  / 0.985, 1.00, 1.01, 0.977, 1.07, 1.31, 1.19, 1.09 /

c                                                             CTEQ1D
      DATA 
     >  Isch(4), Iord(4), Nqrk(4), Alm(4) /  2,  2,  6, .164  / 
     >  (Vlm(I,4), I=4,6) / .247,    .164,    .064  /
     >  (Qms(I,4), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(4), Qmn(4), Qmx(4) /  1.E-5,  2.00,  1.E3  /
     >  Mex(4), Mpn(4), Nexp(4) /  5, 3, 8 /
     >  (ExpN(I, 4), I=1,8 )
     >  / 0.983, 1.00, 1.01, 0.975, 0.964, 1.23, 1.00, 1.12 /
c                                                             CTEQ1L
      DATA 
     >  Isch(5), Iord(5), Nqrk(5), Alm(5) /  1,  1,  6, .125  / 
     >  (Vlm(I,5), I=4,6) / .168,    .125,     .063   /
     >  (Qms(I,5), I=4,6) / 1.50,   5.00,  180.0 /
     >  Xmn(5), Qmn(5), Qmx(5) /  1.E-5,  2.00,  1.E3  /
     >  Mex(5), Mpn(5), Nexp(5) /  5, 3, 8  /
     >  (ExpN(I, 5), I=1,8 )
     >  / 0.982, 1.01, 1.00, 0.972, 0.840, 0.959, 0.930, 0.861 /

      Data ist, lp, qsto, Aln2 / 0, -10, 1.2345, 0.6931 /

      X  = XX
      if(iset.eq.ist.and.iprtn.eq.lp.and.qsto.eq.qq) goto 100

      Irt = 0
      ip = abs(iprtn)
      If (Ip.GE.5.and.QQ.LE.Qms(Ip, Iset)) Then
         CtqPdf = 0.0
         Return
      Endif

c if heavy parton, different logarithmic scale

      if(ip.ge.5) then
         Qi = qms(ip,iset)
      else
         Qi   = Qmn (Iset)
      endif

      Alam = Alm (Iset)
      sta = log(qq/alam)
      stb = log(qi/alam)

c      SBL = LOG(QQ/Alam) / LOG(Qi/Alam)
      sbl = sta/stb
      SB = LOG (SBL)
      SB2 = SB*SB
      SB3 = SB2*SB

      iflv = 3 - iprtn

      Goto (1, 2, 3, 4, 5), Iset

 1    Goto(11,12,13,14,15,16,17,18,19)Iflv    
c   ifl =     2
 11   A0=0.3636E+01*(1.0 + 0.3122E+00*SB+0.1396E+00*SB2+0.4251E+00*SB3)
      A1=0.6930E+00-.2574E-01*SB+0.1047E+00*SB2-.2794E-01*SB3
      A2=0.3195E+01+0.4045E+00*SB-.3737E+00*SB2-.1677E+00*SB3
      A3=0.1009E+00*(1.0 -.1784E+01*SB+0.6263E+00*SB2+0.7337E-01*SB3)
     $  -1.0
      A4=0.2910E+00-.2793E+00*SB+0.6155E-01*SB2+0.5150E-02*SB3
      A5=0.0000E+00+0.3185E+00*SB+0.1953E+00*SB2+0.4184E-01*SB3
      goto 100
c   ifl =     1
 12   A0=0.2851E+00*(1.0 + 0.3617E+00*SB-.4526E+00*SB2+0.5787E-01*SB3)
      A1=0.2690E+00+0.1104E-01*SB+0.1888E-01*SB2-.1031E-01*SB3
      A2=0.3766E+01+0.7850E+00*SB-.3053E+00*SB2+0.1822E+00*SB3
      A3=0.2865E+02*(1.0 -.9774E+00*SB+0.5958E+00*SB2-.1234E+00*SB3)
     $  -1.0
      A4=0.8230E+00-.3612E+00*SB+0.5520E-01*SB2+0.1571E-01*SB3
      A5=0.0000E+00+0.2145E-01*SB+0.2289E+00*SB2-.4947E-01*SB3
      goto 100
c   ifl =     0
 13   A0=0.2716E+01*(1.0 -.2092E+01*SB+0.1500E+01*SB2-.3703E+00*SB3)
      A1=-.3100E-01-.7963E+00*SB+0.1129E+01*SB2-.4191E+00*SB3
      A2=0.8015E+01+0.1168E+01*SB-.1625E+01*SB2-.1130E+01*SB3
      A3=0.4813E+02*(1.0 -.4951E+00*SB-.8715E+00*SB2+0.5893E+00*SB3)
     $  -1.0
      A4=0.2773E+01-.6329E+00*SB-.1048E+01*SB2+0.1418E+00*SB3
      A5=0.0000E+00+0.5048E+00*SB+0.2390E+01*SB2-.4159E+00*SB3
      goto 100
c   ifl =    -1
 14   A0=0.3085E+00*(1.0 + 0.9422E+00*SB-.2606E+01*SB2+0.1364E+01*SB3)
      A1=0.5000E-02-.6433E+00*SB+0.4980E+00*SB2-.1780E+00*SB3
      A2=0.7490E+01+0.9112E+00*SB-.2047E+01*SB2+0.1456E+01*SB3
      A3=0.1145E-01*(1.0 + 0.4610E+01*SB+0.1699E+01*SB2+0.1296E+00*SB3)
     $  -1.0
      A4=0.6030E+00-.8081E+00*SB+0.9410E+00*SB2-.4458E+00*SB3
      A5=0.0000E+00-.1736E+01*SB+0.2863E+01*SB2-.1268E+01*SB3
      goto 100
c   ifl =    -2
 15   A0=0.1324E+00*(1.0 -.1050E+01*SB+0.4844E+00*SB2-.1043E+00*SB3)
      A1=-.1580E+00+0.1672E+00*SB-.4100E+00*SB2+0.1793E+00*SB3
      A2=0.8559E+01-.7351E-01*SB+0.5898E+00*SB2-.2655E+00*SB3
      A3=0.2378E+02*(1.0 -.1108E+00*SB-.1646E-01*SB2+0.1129E-01*SB3)
     $  -1.0
      A4=0.1477E+01+0.3312E-01*SB-.2191E+00*SB2+0.9588E-01*SB3
      A5=0.0000E+00+0.1850E+01*SB-.1481E+01*SB2+0.6222E+00*SB3
      goto 100
c   ifl =    -3
 16   A0=0.3208E+00*(1.0 -.4755E+00*SB-.4003E+00*SB2+0.2300E+00*SB3)
      A1=-.3200E-01-.3357E+00*SB+0.3222E-01*SB2+0.5011E-01*SB3
      A2=0.1164E+02+0.1048E+01*SB-.1097E+01*SB2-.4431E+00*SB3
      A3=0.5065E+02*(1.0 + 0.2484E+00*SB-.9235E+00*SB2+0.1935E+00*SB3)
     $  -1.0
      A4=0.3300E+01-.6785E+00*SB+0.5337E+00*SB2-.4035E+00*SB3
      A5=0.0000E+00-.2496E+00*SB+0.3903E+00*SB2+0.1392E+00*SB3
      goto 100
c   ifl =    -4
 17   A0=0.7967E-06*(1.0 + 0.1587E+01*SB+0.1812E+02*SB2-.1333E+02*SB3)
     $ *sqrt(sta - stb)
      A1=0.1096E+01-.1236E+01*SB+0.1014E+02*SB2+0.1940E+01*SB3
      A2=0.4366E+00+0.1197E+02*SB-.5471E+00*SB2-.5427E+01*SB3
      A3=0.4650E+03*(1.0 + 0.1310E+02*SB-.1918E+02*SB2+0.6791E+01*SB3)
     $  -1.0
      A4=-.8486E+00+0.7457E+00*SB-.1083E+02*SB2-.1210E+01*SB3
      A5=0.3494E+01-.3511E+01*SB-.1766E+01*SB2+0.3442E+01*SB3
      goto 100
c   ifl =    -5
 18   A0=0.1713E-03*(1.0 + 0.2562E+02*SB-.2988E+02*SB2+0.4798E+01*SB3)
     $ *sqrt(sta - stb)
      A1=-.5276E-01+0.4105E+00*SB-.1079E+01*SB2+0.6278E+00*SB3
      A2=0.4515E+01+0.8369E+01*SB-.1192E+02*SB2+0.3403E+01*SB3
      A3=0.1756E+01*(1.0 + 0.1325E+02*SB-.2997E+02*SB2+0.1758E+02*SB3)
     $  -1.0
      A4=0.3557E-01+0.4159E+01*SB-.6947E+01*SB2+0.2982E+01*SB3
      A5=0.2551E+01+0.2168E+01*SB-.5119E+01*SB2+0.3739E+01*SB3
      goto 100
c   ifl =    -6
 19   A0=0.7510E-04*(1.0 + 0.2836E+02*SB-.3000E+02*SB2-.2979E+02*SB3)
     $ *sqrt(sta - stb)
      A1=-.1855E+00+0.4543E+00*SB-.1448E+01*SB2+0.2009E-01*SB3
      A2=0.6775E+01-.4210E+01*SB-.1221E+01*SB2+0.1199E+02*SB3
      A3=0.1070E+01*(1.0 + 0.8356E+01*SB-.2992E+02*SB2+0.2433E+02*SB3)
     $  -1.0
      A4=-.4601E-01+0.4248E+01*SB-.1736E+01*SB2+0.1187E+02*SB3
      A5=0.2771E+01+0.1382E+01*SB-.4797E+01*SB2+0.1273E+01*SB3
      goto 100


 2    Goto(21,22,23,24,25,26,27,28,29)Iflv    
c                                                             CTEQ1MS
c   ifl =     2
 21   A0=0.1828E+01*(1.0 -.8698E+00*SB+0.2906E+00*SB2-.2003E-01*SB3)
      A1=0.6060E+00+0.8595E-01*SB-.4934E-01*SB2+0.2221E-01*SB3
      A2=0.3454E+01-.3115E+00*SB+0.1321E+01*SB2-.3490E+00*SB3
      A3=0.2616E+00*(1.0 -.1670E+01*SB+0.2333E+01*SB2+0.7730E-01*SB3)
     $  -1.0
      A4=0.8920E+00-.8500E-02*SB+0.4960E+00*SB2-.4045E-01*SB3
      A5=0.0000E+00+0.1091E+01*SB-.1613E+00*SB2+0.3773E-01*SB3
      goto 100
c   ifl =     1
 22   A0=0.2885E+00*(1.0 + 0.3388E+00*SB-.4550E+00*SB2+0.6005E-01*SB3)
      A1=0.2730E+00+0.1198E-01*SB+0.1880E-01*SB2-.1077E-01*SB3
      A2=0.3736E+01+0.7687E+00*SB-.2731E+00*SB2+0.1638E+00*SB3
      A3=0.2741E+02*(1.0 -.9585E+00*SB+0.5925E+00*SB2-.1239E+00*SB3)
     $  -1.0
      A4=0.8040E+00-.3546E+00*SB+0.6123E-01*SB2+0.1086E-01*SB3
      A5=0.0000E+00+0.4277E-01*SB+0.2187E+00*SB2-.4646E-01*SB3
      goto 100
c   ifl =     0
 23   A0=0.8416E-01*(1.0 -.1996E+01*SB+0.1903E+01*SB2-.6722E+00*SB3)
      A1=-.4790E+00-.5459E+00*SB+0.1638E+01*SB2-.4342E+00*SB3
      A2=0.5071E+01+0.1470E+01*SB-.2401E+01*SB2+0.1273E+01*SB3
      A3=0.2847E+02*(1.0 + 0.1124E+00*SB-.1338E+01*SB2+0.7115E+00*SB3)
     $  -1.0
      A4=0.4990E+00-.7208E+00*SB+0.3333E-03*SB2-.2354E+00*SB3
      A5=0.0000E+00-.4480E+00*SB+0.3720E+01*SB2-.1838E+01*SB3
      goto 100
c   ifl =    -1
 24   A0=0.4378E+00*(1.0 -.1244E+01*SB+0.3278E+01*SB2-.2098E+01*SB3)
      A1=0.3500E-01-.1298E+01*SB+0.1229E+01*SB2-.3665E+00*SB3
      A2=0.6781E+01+0.4078E+01*SB-.9711E+00*SB2-.1536E+01*SB3
      A3=0.1527E-03*(1.0 + 0.1430E+02*SB+0.3000E+02*SB2+0.2771E+02*SB3)
     $  -1.0
      A4=0.3060E+00+0.1011E+01*SB-.2045E+01*SB2+0.9422E+00*SB3
      A5=0.0000E+00-.3205E+01*SB+0.2683E+01*SB2-.1746E+00*SB3
      goto 100
c   ifl =    -2
 25   A0=0.7413E-01*(1.0 + 0.1291E+01*SB-.2667E+01*SB2+0.1076E+01*SB3)
      A1=-.2730E+00-.1206E+00*SB+0.1828E+00*SB2-.1001E+00*SB3
      A2=0.7719E+01+0.1537E+01*SB-.6410E+00*SB2-.3920E-01*SB3
      A3=0.1799E+02*(1.0 -.1334E+01*SB+0.1916E+01*SB2-.8878E+00*SB3)
     $  -1.0
      A4=0.1167E+01-.9176E-01*SB+0.5132E+00*SB2-.3460E+00*SB3
      A5=0.0000E+00-.5023E+00*SB+0.1951E+01*SB2-.8427E+00*SB3
      goto 100
c   ifl =    -3
 26   A0=0.6551E+00*(1.0 -.5968E-01*SB+0.5621E-02*SB2-.2074E+00*SB3)
      A1=0.2800E-01-.1138E+01*SB+0.1178E+01*SB2-.4425E+00*SB3
      A2=0.7553E+01+0.3996E+01*SB-.4448E+01*SB2+0.1673E+01*SB3
      A3=0.9264E-01*(1.0 -.1760E+01*SB+0.1634E+01*SB2-.4067E+00*SB3)
     $  -1.0
      A4=0.1970E+00+0.5256E+00*SB-.9775E+00*SB2+0.4488E+00*SB3
      A5=0.0000E+00-.3668E+01*SB+0.4757E+01*SB2-.1717E+01*SB3
      goto 100
c   ifl =    -4
 27   A0=0.1486E-03*(1.0 + 0.2107E+01*SB-.1056E+02*SB2+0.1403E+02*SB3)
     $ * sqrt(sta - stb)
      A1=0.2115E+00-.1702E+01*SB+0.2571E+01*SB2-.1177E+01*SB3
      A2=0.3533E+01+0.1367E+01*SB-.3397E+01*SB2+0.6260E+01*SB3
      A3=0.1096E+02*(1.0 + 0.9213E+01*SB-.2020E+02*SB2+0.1084E+02*SB3)
     $  -1.0
      A4=0.7041E+00-.7236E+00*SB+0.2766E-01*SB2+0.7352E+00*SB3
      A5=0.3904E+01-.4398E+01*SB+0.7056E+01*SB2-.3722E+01*SB3
      goto 100
c   ifl =    -5
 28   A0=0.1201E-03*(1.0 + 0.5408E+01*SB-.1489E+02*SB2+0.1667E+02*SB3)
     $  * sqrt(sta - stb)
      A1=0.1420E-01-.1525E+01*SB+0.2408E+01*SB2-.1154E+01*SB3
      A2=0.4254E+01+0.2836E+01*SB-.6018E+00*SB2+0.4133E+00*SB3
      A3=0.5696E+01*(1.0 + 0.9451E+01*SB-.2029E+02*SB2+0.1033E+02*SB3)
     $  -1.0
      A4=0.4775E+00-.6695E+00*SB+0.2747E+00*SB2-.1051E+00*SB3
      A5=0.3330E+01-.5133E+01*SB+0.6921E+01*SB2-.3283E+01*SB3
      goto 100
c   ifl =    -6
 29   A0=0.7697E-04*(1.0 + 0.2801E+02*SB-.1901E+02*SB2-.2880E+02*SB3)
     $ *sqrt(sta - stb)
      A1=-.2249E+00+0.4432E+00*SB-.1454E+01*SB2+0.3509E-01*SB3
      A2=0.6642E+01-.2702E+01*SB+0.8229E+01*SB2+0.8243E+01*SB3
      A3=0.1146E+01*(1.0 + 0.8104E+01*SB-.2998E+02*SB2+0.2812E+02*SB3)
     $  -1.0
      A4=-.6421E-01+0.4246E+01*SB-.2908E+01*SB2+0.9686E-02*SB3
      A5=0.2606E+01+0.1261E+01*SB-.4933E+01*SB2+0.3476E+00*SB3
      goto 100


 3    Goto(31,32,33,34,35,36,37,38,39)Iflv    
c                                                             CTEQ1ML
c   ifl =     2
 31   A0=0.3777E+01*(1.0 + 0.6986E+00*SB-.20655E+01*SB2+.10334E+01*SB3)
      A1=0.7100E+00+.2880E-01*SB-.7930E-01*SB2+0.5600E-01*SB3
      A2=0.3259E+01+0.1508E+01*SB-.3932E+01*SB2+0.20613E+01*SB3
      A3=0.1304E+00*(1.0 -.2016E+00*SB-.30015E+01*SB2+0.19118E+01*SB3)
     $     -1.0
      A4=0.2890E+00-0.4311E+00*SB+0.7387E+00*SB2-.3697E+00*SB3
      A5=0.0000E+00+0.4320E+00*SB+0.2449E+00*SB2-0.6670E-01*SB3
      goto 100
c   ifl =     1
 32   A0=0.2780E+00*(1.0 + 0.4355E+00*SB-0.4584E+00*SB2+0.4390E-01*SB3)
      A1=0.2760E+00+0.1420E-01*SB+0.1480E-01*SB2-.9800E-02*SB3
      A2=0.3710E+01+0.8250E+00*SB-.3581E+00*SB2+0.1978E+00*SB3
      A3=0.2928E+02*(1.0 -.10154E+01*SB+0.6037E+00*SB2-.1175E+00*SB3)
     $     -1.0
      A4=0.8070E+00-.3575E+00*SB+0.4920E-01*SB2+0.1584E-01*SB3
      A5=0.0000E+00+0.1860E-01*SB+0.2080E+00*SB2-.450E-01*SB3
      goto 100
c   ifl =     0
 33   A0=0.2924E+01*(1.0 -.18916E+01*SB+0.1191E+01*SB2-.2492E+00*SB3)
      A1=0.0000E+00-.9167E+00*SB+0.11147E+01*SB2-.3329E+00*SB3
      A2=0.8529E+01+0.7080E+00*SB-.11345E+01*SB2-.10563E+01*SB3
      A3=0.1420E+03*(1.0 -.15346E+01*SB+0.7261E+00*SB2-.5730E-01*SB3)
     $     -1.0
      A4=0.3396E+01-.11541E+01*SB-.8834E+00*SB2+0.2430E+00*SB3
      A5=0.0000E+00+0.1645E+00*SB+0.19041E+01*SB2+0.1474E+00*SB3
      goto 100
c   ifl =    -1
 34   A0=0.3471E+00*(1.0- 0.1753E+00*SB-.9189E+00*SB2+0.6211E+00*SB3)
      A1=0.1900E-01-.4579E+00*SB+0.2112E+00*SB2-.6180E-01*SB3
      A2=0.7301E+01-.17308E+01*SB+.13666E+01*SB2-.6400E-02*SB3
      A3=0.1853E-04*(1.0 -.18260E+02*SB-.2872E+02*SB2-.23456E+02*SB3)
     $     -1.0
      A4=0.4400E+00-.4672E+00*SB+0.6532E+00*SB2-.3222E+00*SB3
      A5=0.0000E+00-.4679E+00*SB+0.10741E+01*SB2-.5663E+00*SB3
      goto 100
c   ifl =    -2
 35   A0=0.1702E+00*(1.0 -.1041E+01*SB+0.4064E+00*SB2-.5888E-01*SB3)
      A1=-.9300E-01-.4742E-01*SB-.1959E+00*SB2+0.1039E+00*SB3
      A2=0.9119E+01-.7331E-01*SB+0.3506E+00*SB2-.2081E+00*SB3
      A3=0.2981E+02*(1.0 -.1912E+00*SB-.8947E-02*SB2+0.8805E-02*SB3)
     $     -1.0
      A4=0.1668E+01-.6678E-02*SB-.2894E+00*SB2+0.1221E+00*SB3
      A5=0.0000E+00+0.1245E+01*SB-.7843E+00*SB2+0.3724E+00*SB3
      goto 100
c   ifl =    -3
 36   A0=0.3910E+00*(1.0 -.1103E+01*SB+0.5383E+00*SB2-.1083E+00*SB3)
      A1=-.1400E-01-.2471E+00*SB-.8042E-01*SB2+0.7193E-01*SB3
      A2=0.9812E+01-.4860E+01*SB+0.5958E+01*SB2-.2342E+01*SB3
      A3=0.3749E+00*(1.0 -.3569E+01*SB+0.5456E+01*SB2-.2344E+01*SB3)
     $     -1.0
      A4=0.4940E+00+0.2772E+00*SB-.2732E+00*SB2+0.6466E-01*SB3
      A5=0.0000E+00+0.3927E+00*SB-.3216E+00*SB2+0.2164E+00*SB3
      goto 100
c   ifl =    -4
 37   A0=0.3815E-02*(1.0 + 0.2039E+02*SB-.2834E+02*SB2+0.1070E+02*SB3)
     $ * sqrt(sta - stb)
      A1=-.2789E-01-.7345E-03*SB-.3251E+00*SB2+0.1946E+00*SB3
      A2=0.3223E+01-.4268E+00*SB+0.4387E+01*SB2-.2401E+01*SB3
      A3=0.3338E-01*(1.0 -.1163E+02*SB+0.2995E+02*SB2-.1471E+02*SB3)
     $     -1.0
      A4=0.3646E+00-.5767E+00*SB+0.6088E+00*SB2-.2514E+00*SB3
      A5=0.1200E+01+0.2178E+00*SB-.4230E+00*SB2+0.4739E+00*SB3
      goto 100
c   ifl =    -5
 38   A0=0.1666E-02*(1.0 + 0.9518E+01*SB-.4715E+01*SB2-.1060E+01*SB3)
     $ * sqrt(sta - stb)
      A1=-.1231E+00+0.1656E+00*SB-.5219E+00*SB2+0.2750E+00*SB3
      A2=0.3693E+01+0.4922E+01*SB-.1200E+02*SB2+0.7929E+01*SB3
      A3=0.1778E+00*(1.0 + 0.3036E+01*SB-.1184E+02*SB2+0.7940E+01*SB3)
     $     -1.0
      A4=0.5353E+00-.1401E+01*SB+0.1970E+01*SB2-.9405E+00*SB3
      A5=0.1590E+01+0.1025E+01*SB-.2318E+01*SB2+0.1380E+01*SB3
      goto 100
c   ifl =    -6
 39   A0=0.4319E-03*(1.0 + 0.1100E+02*SB-.9520E+00*SB2+0.1434E+02*SB3)
     $ * sqrt(sta - stb)
      A1=-.2512E+00+0.3554E+00*SB-.4120E+00*SB2+0.1328E+00*SB3
      A2=0.4764E+01-.3513E+00*SB+0.1199E+02*SB2-.8290E+01*SB3
      A3=0.8458E-01*(1.0 + 0.2618E+01*SB+0.4407E+01*SB2+0.2991E+02*SB3)
     $    -1.0
      A4=0.3991E+00-.1363E+01*SB+0.1526E+01*SB2-.3179E+01*SB3
      A5=0.1981E+01+0.1496E+01*SB-.1501E+01*SB2+0.3880E+01*SB3
      goto 100

 4    Goto(41,42,43,44,45,46,47,48,49)Iflv
c                                                             CTEQ1D
c   ifl =     2
 41   A0=0.1634E+01*(1.0 -.8336E+00*SB+0.1640E+00*SB2+0.1530E+00*SB3)
      A1=0.5790E+00+0.8587E-01*SB-.6087E-01*SB2+0.1361E-01*SB3
      A2=0.2839E+01+0.3720E+00*SB+0.5264E+00*SB2+0.3538E-01*SB3
      A3=0.1095E+00*(1.0 -.4830E+00*SB+0.3708E+01*SB2-.6165E+00*SB3)
     $  -1.0
      A4=0.8010E+00-.1432E+00*SB+0.1442E+01*SB2-.1286E+01*SB3
      A5=0.0000E+00+0.1035E+01*SB-.5910E-01*SB2-.1982E+00*SB3
      goto 100
c   ifl =     1
 42   A0=0.3535E+00*(1.0 + 0.4352E+00*SB-.2095E+00*SB2-.8455E-02*SB3)
      A1=0.2660E+00-.4096E-03*SB+0.1502E-01*SB2-.1163E-01*SB3
      A2=0.3514E+01+0.8219E+00*SB-.2330E+00*SB2+0.1055E+00*SB3
      A3=0.2200E+02*(1.0 -.9716E+00*SB+0.4552E+00*SB2-.8202E-01*SB3)
     $  -1.0
      A4=0.9000E+00-.3207E+00*SB-.4808E-01*SB2+0.3492E-01*SB3
      A5=0.0000E+00-.6273E-01*SB+0.1497E+00*SB2-.5683E-01*SB3
      goto 100
c   ifl =     0
 43   A0=0.2743E+01*(1.0 -.2027E+01*SB+0.1517E+01*SB2-.4145E+00*SB3)
      A1=0.7000E-02-.9431E+00*SB+0.1231E+01*SB2-.4834E+00*SB3
      A2=0.8200E+01+0.1827E+01*SB-.3453E+01*SB2+0.6763E+00*SB3
      A3=0.4975E+02*(1.0 -.2329E+00*SB-.1245E+01*SB2+0.7194E+00*SB3)
     $  -1.0
      A4=0.2387E+01-.4077E+00*SB-.5542E+00*SB2-.9677E-02*SB3
      A5=0.0000E+00+0.2702E+00*SB+0.2389E+01*SB2-.8274E+00*SB3
      goto 100
c   ifl =    -1
 44   A0=0.2015E+00*(1.0 -.2133E+00*SB-.6770E+00*SB2+0.5011E+00*SB3)
      A1=-.7700E-01-.7104E-01*SB-.3720E+00*SB2+0.2159E+00*SB3
      A2=0.8008E+01-.2049E+01*SB+0.1800E+01*SB2-.4660E+00*SB3
      A3=0.2923E-05*(1.0 + 0.2327E+02*SB+0.1500E+02*SB2+0.2633E+02*SB3)
     $  -1.0
      A4=0.9020E+00-.9191E+00*SB+0.1104E+01*SB2-.5863E+00*SB3
      A5=0.0000E+00+0.5840E+00*SB-.8720E+00*SB2+0.4234E+00*SB3
      goto 100
c   ifl =    -2
 45   A0=0.9117E-01*(1.0 -.4089E+00*SB-.4361E+00*SB2+0.2512E+00*SB3)
      A1=-.2370E+00+0.2492E+00*SB-.3267E+00*SB2+0.1055E+00*SB3
      A2=0.8447E+01+0.6009E+00*SB+0.1003E+01*SB2-.1287E+01*SB3
      A3=0.3106E+02*(1.0 -.3901E-01*SB+0.1443E+00*SB2-.3433E+00*SB3)
     $  -1.0
      A4=0.1629E+01+0.7855E-01*SB-.1573E+00*SB2-.8595E-01*SB3
      A5=0.0000E+00+0.1558E+01*SB-.6295E+00*SB2+0.1847E+00*SB3
      goto 100
c   ifl =    -3
 46   A0=0.3997E+00*(1.0 -.1046E+01*SB+0.6194E+00*SB2-.1342E+00*SB3)
      A1=0.2000E-02-.2544E+00*SB-.1958E+00*SB2+0.1458E+00*SB3
      A2=0.9613E+01-.3919E+01*SB+0.9573E+01*SB2-.5623E+01*SB3
      A3=0.3620E+00*(1.0 -.1858E+01*SB+0.8312E+01*SB2-.5900E+01*SB3)
     $  -1.0
      A4=0.3840E+00+0.3572E+00*SB-.1191E+01*SB2+0.7310E+00*SB3
      A5=0.0000E+00+0.3351E+00*SB-.7709E+00*SB2+0.4296E+00*SB3
      goto 100
c   ifl =    -4
 47   A0=0.2156E-03*(1.0 + 0.2879E+02*SB-.2310E+02*SB2+0.9812E+01*SB3)
     $ * sqrt(sta - stb)
      A1=0.9086E-01-.1250E+00*SB-.7373E-01*SB2-.2201E-01*SB3
      A2=0.3588E+01+0.4518E+01*SB-.8930E-01*SB2+0.9163E-02*SB3
      A3=0.5216E+01*(1.0 + 0.5912E+00*SB-.4111E+00*SB2+0.7330E+00*SB3)
     $  -1.0
      A4=0.3145E+00+0.1233E+01*SB-.7478E+00*SB2+0.4657E+00*SB3
      A5=0.2723E+01-.4110E+00*SB+0.4868E-01*SB2-.3075E+00*SB3
      goto 100
c   ifl =    -5
 48   A0=0.7476E-03*(1.0 + 0.1454E+02*SB-.2509E+02*SB2+0.1184E+02*SB3)
     $ * sqrt(sta - stb)
      A1=-.1955E-01-.1712E+00*SB-.1686E+00*SB2+0.2339E+00*SB3
      A2=0.4616E+01-.6859E+00*SB-.3959E+01*SB2+0.5530E+01*SB3
      A3=0.9881E+01*(1.0 -.1239E+02*SB+0.2721E+02*SB2-.1850E+02*SB3)
     $  -1.0
      A4=0.1200E+02-.1133E+02*SB+0.8138E+01*SB2+0.1199E+02*SB3
      A5=0.2226E+01-.5738E+00*SB+0.5239E+00*SB2+0.3825E+00*SB3
      goto 100
c   ifl =    -6
 49   A0=0.8392E-06*(1.0 + 0.1844E+02*SB-.1110E+02*SB2-.2504E+02*SB3)
     $ * sqrt(sta - stb)
      A1=0.2127E+00-.5602E+00*SB+0.4777E+01*SB2-.1014E+02*SB3
      A2=0.1229E+01+0.7495E+01*SB-.5024E+01*SB2-.1200E+02*SB3
      A3=0.2868E+02*(1.0 + 0.7634E+01*SB-.2916E+02*SB2+0.2953E+02*SB3)
     $  -1.0
      A4=0.5970E+00+0.1138E+01*SB-.1439E+01*SB2-.1966E+01*SB3
      A5=0.6429E+01-.6673E+00*SB+0.7008E+01*SB2-.1157E+02*SB3
      goto 100

 5    Goto(51,52,53,54,55,56,57,58,59)Iflv
c                                                             CTEQ1L
c   ifl =     2
 51   A0=  1.791*(1.0 -0.449*SB-0.445*SB2+  0.401*SB3)
      A1=  0.608+  0.069*SB+  0.005*SB2-0.037*SB3
      A2=  3.470-0.375*SB+  2.267*SB2-1.261*SB3
      A3=  0.315*(1.0 -2.628*SB+  6.481*SB2-3.834*SB3)-1.0
      A4=  1.007-0.732*SB+  1.490*SB2-0.966*SB3
      A5=  0.000+  0.741*SB+  0.563*SB2-0.525*SB3
      goto 100
c   ifl =     1
 52   A0=  0.513*(1.0 +   0.032*SB-0.120*SB2+  0.013*SB3)
      A1=  0.276+  0.052*SB+  0.000*SB2-0.006*SB3
      A2=  3.579+  0.763*SB-0.135*SB2+  0.083*SB3
      A3= 17.993*(1.0 -0.725*SB+  0.241*SB2-0.020*SB3)-1.0
      A4=  1.120-0.357*SB+  0.008*SB2+  0.028*SB3
      A5=  0.000+  0.311*SB+  0.029*SB2-0.010*SB3
      goto 100
c   ifl =     0
 53   A0=  2.710*(1.0 -1.773*SB+  0.970*SB2-0.149*SB3)
      A1= -0.010-1.636*SB+  2.087*SB2-0.637*SB3
      A2=  7.174+  2.102*SB-2.209*SB2-0.420*SB3
      A3= 29.904*(1.0 -0.756*SB-0.506*SB2+  0.605*SB3)-1.0
      A4=  2.572-0.437*SB-0.968*SB2+  0.243*SB3
      A5=  0.000-1.776*SB+  4.266*SB2-0.335*SB3
      goto 100
c   ifl =    -1
 54   A0=  0.278*(1.0 - 1.022*SB+  0.6457*SB2-0.1824*SB3)
      A1=  0.0862*SB-0.8657*SB2+  0.4185*SB3
      A2= 11.000-1.2809*SB+ 1.2516*SB2+0.061*SB3
      A3= 37.338*(1.0 - 0.9404*SB+  0.2517*SB2+0.1364*SB3)-1.0
      A4=  1.960-  0.3385*SB-0.3422*SB2+0.3653*SB3
      A5=  0.000+1.424*SB-2.7503*SB2+  1.2226*SB3
      goto 100
c   ifl =    -2
 55   A0=  0.154*(1.0 -0.659*SB+  0.005*SB2+  0.061*SB3)
      A1= -0.128+  0.279*SB-0.786*SB2+  0.363*SB3
      A2=  8.649+  0.071*SB+  0.351*SB2-0.051*SB3
      A3= 43.685*(1.0 -0.603*SB+  0.037*SB2+  0.134*SB3)-1.0
      A4=  2.238-0.338*SB-0.199*SB2+  0.157*SB3
      A5=  0.000+  1.681*SB-2.068*SB2+  0.975*SB3
      goto 100
c   ifl =    -3
 56   A0=  0.372*(1.0 -1.939*SB+  1.504*SB2-0.440*SB3)
      A1=  0.009+  0.610*SB-1.387*SB2+  0.579*SB3
      A2= 10.273-4.833*SB+  6.583*SB2-2.633*SB3
      A3=  0.160*(1.0 +  10.325*SB-2.027*SB2+  1.571*SB3)-1.0
      A4=  0.819-1.660*SB+  1.845*SB2-0.829*SB3
      A5=  0.000+  3.558*SB-3.940*SB2+  1.302*SB3
      goto 100
c   ifl =    -4
 57   A0=  (7.5242E-5)*(1.0+22.0905*SB+7.1209*SB2-8.303*SB3)*
     $     sqrt(sta - stb)
      A1=  0.125-0.3027*SB+0.1564*SB2-0.091*SB3
      A2=  2.0388+1.2161*SB+11.5296*SB2-8.0659*SB3
      A3=  14.849*(1.0 -2.556*SB+3.5268*SB2-1.6353*SB3)-1.0
      A4=  0.3061-0.0901*SB+0.953*SB2-0.4871*SB3
      A5=  2.7352+0.1811*SB-0.5167*SB2+0.0543*SB3
      goto 100
c   ifl =    -5
 58   A0=  (3.751E-4)*(1.0 + 21.5993*SB+3.1379*SB2-18.8328*SB3)*
     $     sqrt(sta - stb)
      A1= -0.0256-0.7717*SB+ 1.1499*SB2-0.5037*SB3
      A2=  4.9241+4.0107*SB-4.7012*SB2+0.1097*SB3
      A3=  2.842*(1.0 -2.2184*SB+  2.0293*SB2-0.6907*SB3)-1.0
      A4=  -0.1352+ 0.8753*SB-1.2626*SB2+  0.667*SB3
      A5=  1.5627-0.4917*SB+ 1.5927*SB2-0.351*SB3
      goto 100
c   ifl =    -6
 59   A0=(2.725E-4)*(1.0 +  18.8497*SB-26.5797*SB2-29.0774*SB3)*
     $     sqrt(sta - stb)
      A1= -0.2204-1.0048*SB+0.9415*SB2-0.4274*SB3
      A2=  11.034-9.8362*SB-11.1034*SB2-9.1977*SB3
      A3=  2.084*(1.0 -2.881*SB+1.2778*SB2-2.9328*SB3)-1.0
      A4= -0.0872+  0.200*SB-1.6187*SB2-1.6058*SB3
      A5=  0.8684+4.7047*SB-1.4614*SB2-5.2309*SB3
      goto 100

 100  CtqPdf = A0*(x**A1)*((1.-x)**A2)*(1.+A3*(x**A4))
     $     *(log(1.+1./x))**A5

      if(ctqpdf.lt.0.0) ctqpdf = 0.0

      Ist = Iset

      Lp  = Iprtn
      Qsto = QQ

      Return
C                                  -----------------------
      ENTRY WLAMBD (ISET, IORDER)

      IORDER = IORD (ISET)
      WLAMBD = ALM  (ISET)

      RETURN
C                                  -----------------------
      Entry PrCtq1 
     >        (Iset, Iordr, Ischeme, MxFlv,
     >         Alam4, Alam5, Alam6, Amas4, Amas5, Amas6,
     >         Xmin, Qini, Qmax, ExpNor)

C                           Return QCD parameters and Fitting parameters
C                           associated with parton distribution set Iset.
C    Iord    : Order Of Fit
C    Ischeme : (0, 1, 2)  for  (LO, MS-bar-NLO, DIS-NLO) resp.
C    MxFlv   : Maximum number of flavors included
C    Alam_i  : i = 4,5,6  Effective lambda for i-flavors 

C    Amas_i  : i = 4,5,6  Mass parameter for flavor i
C    Xmin, Qini, Qmax : self explanary
C    ExpNor(I) : Normalization factor for the experimental data set used in
C                obtaining the best global fit for parton distributions Iset:
C     I = 1,     2,      3,     4,     5,     6,     7,     8
C      BCDMS   NMC90  NMC280  CCFR   E605    WA70   E706   UA6

      Iordr  = Iord (Iset)
      Ischeme= Isch (Iset)
      MxFlv  = Nqrk (Iset)

      Alam4  = Vlm(4,Iset)
      Alam5  = Vlm(5,Iset)
      Alam6  = Vlm(6,Iset)

      Amas4  = Qms(4,Iset)
      Amas5  = Qms(5,Iset)
      Amas6  = Qms(6,Iset)

      Xmin   = Xmn  (Iset)
      Qini   = Qmn  (Iset)
      Qmax   = Qmx  (Iset)

      Do 101 Iexp = 1, Nexp(Iset)
         ExpNor(Iexp) = ExpN(Iexp, Iset)
  101 Continue

      Return
C                         *************************
      END


