MODULE KPP_ROOT_Main
  ! A driver option for calling KPP as a module, from an external program
  ! by Sampo Smolander, 2008

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: KPP_SetUp, KPP_Proceed, PHOTOLYSIS_KPP

  KPP_REAL :: T, DVAL(NSPEC)
  KPP_REAL :: RSTATE(20)
  INTEGER :: ICNTRL(20) ! Added by Sampo Smolander
  INTEGER :: i

CONTAINS
  
!~~~> Initialization 

  SUBROUTINE KPP_SetUp()
    STEPMIN = 0.0d0
    STEPMAX = 0.0d0
    ICNTRL = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
    DO i=1,NVAR
       RTOL(i) = 1.0d-4
       ATOL(i) = 1.0d-3
    END DO
    CALL Initialize()
  END SUBROUTINE KPP_SetUp

!~~~> Time loop

  SUBROUTINE KPP_Proceed(CONS,TSTARTp,TENDp,TEMPp,RHp,PRESSp,O2p,N2p,Mp,H2Op,RES1p,RES2p,Jp)

    IMPLICIT NONE
    REAL(kind=dp), INTENT(INOUT) :: CONS(NSPEC)
    REAL(kind=dp), INTENT(IN) :: TSTARTp,TENDp
    REAL(kind=dp), INTENT(IN) :: TEMPp,RHp,PRESSp,O2p,N2p,Mp,H2Op,RES1p,RES2p
    REAL(kind=dp), INTENT(IN) :: Jp(NPHOT)

    REAL(kind=dp) :: T1,T2

    T1 = TSTARTp
    T2 = TENDp

    TEMP  = TEMPp
    RH    = RHp
    PRESS = PRESSp
    O2    = O2p
    N2    = N2p
    M     = Mp
    H2O   = H2Op
    RES1  = RES1p
    RES2  = RES2p
    J     = Jp

    ! Update_RCONST and INTEGRATE will operate on C (actually on VAR and FIX)
    C = CONS

    CALL Update_RCONST()

    ! In INTEGRATE VAR and FIX are calculated, so we need to put values from C to them
    ! And FIX is meaningful only when NVAR < NSPEC
    VAR(1:NVAR) = C(1:NVAR)
    DO i = NVAR+1, NSPEC
      FIX(i-NVAR) = C(i)
    END DO

    CALL INTEGRATE( TIN = T1, TOUT = T2, RSTATUS_U = RSTATE, ICNTRL_U = ICNTRL )

    ! Set VAR and FIX back to C
    C(1:NVAR) = VAR(1:NVAR)
    DO i = NVAR+1, NSPEC
      C(i) = FIX(i-NVAR)
    END DO

    CONS = C

  END SUBROUTINE KPP_Proceed


  SUBROUTINE PHOTOLYSIS_KPP(Beta, SAFN, SAFN_out, T, Jpr)

    IMPLICIT NONE

    INTEGER ::                      I

    INTEGER, DIMENSION(75), SAVE :: WL !wavelength

    REAL(kind=dp), DIMENSION(75) :: SAFN, SAFN_out ! Actinic flux in two different units

    REAL(kind=dp) :: T,                                                                 & ! temperature in K
         Beta,                                                                          & ! solar zenit angle
         DL,                                                                            & ! Delta lambda (so wl band)
         VAA, VAB, VAC, VAD, h_planck, c_light,                                         & ! Dimension parameter
         Jpr(NPHOT),                                                                    & ! photolysis rate 
         ACS_o3, ACS_hcho, ACS_ch3coch3, ACS_ch3no3, ACS_c2h5no3, ACS_i_c3h7no3           ! some cross sections are more comples

    REAL(kind=dp), DIMENSION(75), SAVE :: ACS_mek, ACS_h2o2, ACSA_o3, ACSB_o3, ACS_no2, & ! ACS_X : absorption cross section
         QY_mek, QY_h2o2, QY_o1d, QY_o3p, QY_no2, QY_hono, ACS_hono, QY_hno3, ACS_hno3, & ! QY_x : quantum yield for photolysis
         ACS_no3, ACSA_hcho, QY_no3_no, QY_no3_no2, ACSB_hcho, QY_hcho_h_hco,           &
         QY_hcho_h2_co, ACS_ch3cho, QY_ch3cho, ACS_c2h5cho, QY_c2h5cho, ACS_ic3h7cho,   &
         QY_ic3h7cho, ACS_nc3h7cho, QY1_nc3h7cho, QY2_nc3h7cho, ACS_macr, QY1_macr,     &
         QY2_macr, ACSA_ch3coch3, ACSB_ch3coch3, ACSC_ch3coch3, ACSD_ch3coch3,          &
         QY_ch3coch3, ACS_mvk, QY1_mvk, QY2_mvk, ACS_glyox, QY1_glyox, QY3_glyox,       &
         QY2_glyox, ACS_mglyox, QY_mglyox, ACS_biacet, QY_biacet, ACS_ch3ooh,           &
         QY_ch3ooh, ACSA_ch3no3, ACSB_ch3no3, QY_ch3no3, ACSA_c2h5no3, ACSB_c2h5no3,    &
         QY_c2h5no3, ACS_n_c3h7no3, QY_n_c3h7no3, QY_i_c3h7no3, ACSA_i_c3h7no3,         &
         ACSB_i_c3h7no3, ACS_t_c4h9no3, QY_t_c4h9no3, ACS_noa, QY1_noa, QY2_noa,        &
         ACS_ho2no2, QY1_ho2no2, QY2_ho2no2, ACS_n2o5, ACSA_n2o5, ACSB_n2o5, QY1_n2o5,  &
         QY2_n2o5, ACS_test

    REAL :: VAAAA, T1, AAAA, BBBB, CCCC, LOOO, QY_o3a, KULMA

    LOGICAL, SAVE :: first_call = .TRUE.

    CHARACTER(*), PARAMETER :: filename1 = '/home/local/pzzhou/Scripts/Fortran/SOSAA/sosa_in'

    ! The following is run only once, only on the first time CHEMISTRY is called
    ! Input of absorption cross spectrum and quantum yields for different molecules: 
    IF (first_call) THEN ! do only once
      first_call = .FALSE.

      OPEN(900,  FILE = ''//filename1//'/General/Photolyse/o3_mcm_v2.dat',          STATUS = 'OLD')
      OPEN(901,  FILE = ''//filename1//'/General/Photolyse/o1d_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(902,  FILE = ''//filename1//'/General/Photolyse/o3p_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(903,  FILE = ''//filename1//'/General/Photolyse/h2o2_mcm_v2.dat',        STATUS = 'OLD')
      OPEN(904,  FILE = ''//filename1//'/General/Photolyse/no2_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(905,  FILE = ''//filename1//'/General/Photolyse/no2_qy_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(906,  FILE = ''//filename1//'/General/Photolyse/no3_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(907,  FILE = ''//filename1//'/General/Photolyse/no3_no_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(908,  FILE = ''//filename1//'/General/Photolyse/no3_no2_mcm_v2.dat',     STATUS = 'OLD')
      OPEN(909,  FILE = ''//filename1//'/General/Photolyse/hono_mcm_v2.dat',        STATUS = 'OLD')
      OPEN(910,  FILE = ''//filename1//'/General/Photolyse/hno3_mcm_v2_298.dat',    STATUS = 'OLD')
      OPEN(911,  FILE = ''//filename1//'/General/Photolyse/hcho_mcm_v2.dat',        STATUS = 'OLD')
      OPEN(912,  FILE = ''//filename1//'/General/Photolyse/hcho_h_hco_mcm_v2.dat',  STATUS = 'OLD')
      OPEN(913,  FILE = ''//filename1//'/General/Photolyse/hcho_h2_co_mcm_v2.dat',  STATUS = 'OLD')
      OPEN(914,  FILE = ''//filename1//'/General/Photolyse/ch3cho_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(915,  FILE = ''//filename1//'/General/Photolyse/ch3cho_qy_mcm_v2.dat',   STATUS = 'OLD')
      OPEN(916,  FILE = ''//filename1//'/General/Photolyse/c2h5cho_mcm_v2.dat',     STATUS = 'OLD')
      OPEN(917,  FILE = ''//filename1//'/General/Photolyse/c2h5cho_qy_mcm_v2.dat',  STATUS = 'OLD')        
      OPEN(918,  FILE = ''//filename1//'/General/Photolyse/nc3h7cho_mcm_v2.dat',    STATUS = 'OLD')
      OPEN(919,  FILE = ''//filename1//'/General/Photolyse/ic3h7cho_mcm_v2.dat',    STATUS = 'OLD')
      OPEN(920,  FILE = ''//filename1//'/General/Photolyse/ic3h7cho_qy_mcm_v2.dat', STATUS = 'OLD')
      OPEN(921,  FILE = ''//filename1//'/General/Photolyse/macr_mcm_v2.dat',        STATUS = 'OLD')
      OPEN(922,  FILE = ''//filename1//'/General/Photolyse/ch3coch3_mcm_v2.dat',    STATUS = 'OLD')
      OPEN(923,  FILE = ''//filename1//'/General/Photolyse/ch3coch3_qy_mcm_v2.dat', STATUS = 'OLD')
      OPEN(924,  FILE = ''//filename1//'/General/Photolyse/mek_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(925,  FILE = ''//filename1//'/General/Photolyse/mvk_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(926,  FILE = ''//filename1//'/General/Photolyse/mvk_qy_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(927,  FILE = ''//filename1//'/General/Photolyse/glyox_mcm_v2.dat',       STATUS = 'OLD')
      OPEN(928,  FILE = ''//filename1//'/General/Photolyse/glyox_qy_mcm_v2.dat',    STATUS = 'OLD')
      OPEN(929,  FILE = ''//filename1//'/General/Photolyse/mglyox_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(930,  FILE = ''//filename1//'/General/Photolyse/mglyox_qy_mcm_v2.dat',   STATUS = 'OLD')
      OPEN(931,  FILE = ''//filename1//'/General/Photolyse/biacet_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(932,  FILE = ''//filename1//'/General/Photolyse/ch3ooh_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(933,  FILE = ''//filename1//'/General/Photolyse/ch3no3_mcm_v2.dat',      STATUS = 'OLD')
      OPEN(934,  FILE = ''//filename1//'/General/Photolyse/c2h5no3_mcm_v2.dat',     STATUS = 'OLD')
      OPEN(935,  FILE = ''//filename1//'/General/Photolyse/n_c3h7no3_mcm_v2.dat',   STATUS = 'OLD')
      OPEN(936,  FILE = ''//filename1//'/General/Photolyse/i_c3h7no3_mcm_v2.dat',   STATUS = 'OLD')
      OPEN(937,  FILE = ''//filename1//'/General/Photolyse/tc4h9no3_mcm_v2.dat',    STATUS = 'OLD')
      OPEN(938,  FILE = ''//filename1//'/General/Photolyse/noa_mcm_v2.dat',         STATUS = 'OLD')
      OPEN(939,  FILE = ''//filename1//'/General/Photolyse/ho2no2_aitkinson.dat',   STATUS = 'OLD')
      OPEN(940,  FILE = ''//filename1//'/General/Photolyse/n2o5_aitkinson.dat',     STATUS = 'OLD')
      OPEN(941,  FILE = ''//filename1//'/General/Photolyse/test_O3.dat',            STATUS = 'OLD')

      DO I=1,75
         READ(900,*)  WL(I), ACSA_o3(I), ACSB_o3(I)
         READ(901,*)  WL(I), QY_o1d(I)
         READ(902,*)  WL(I), QY_o3p(I)
         READ(903,*)  WL(I), ACS_h2o2(I), QY_h2o2(I)
         READ(904,*)  WL(I), ACS_no2(I)
         READ(905,*)  WL(I), QY_no2(I)
         READ(906,*)  WL(I), ACS_no3(I)
         READ(907,*)  WL(I), QY_no3_no(I)
         READ(908,*)  WL(I), QY_no3_no2(I)
         READ(909,*)  WL(I), ACS_hono(I), QY_hono(I)
         READ(910,*)  WL(I), ACS_hno3(I), QY_hno3(I)
         READ(911,*)  WL(I), ACSA_hcho(I), ACSB_hcho(I)
         READ(912,*)  WL(I), QY_hcho_h_hco(I)
         READ(913,*)  WL(I), QY_hcho_h2_co(I)
         READ(914,*)  WL(I), ACS_ch3cho(I)
         READ(915,*)  WL(I), QY_ch3cho(I)
         READ(916,*)  WL(I), ACS_c2h5cho(I)
         READ(917,*)  WL(I), QY_c2h5cho(I)
         READ(918,*)  WL(I), ACS_nc3h7cho(I), QY1_nc3h7cho(I), QY2_nc3h7cho(I)
         READ(919,*)  WL(I), ACS_ic3h7cho(I)
         READ(920,*)  WL(I), QY_ic3h7cho(I)
         READ(921,*)  WL(I), ACS_macr(I), QY1_macr(I), QY2_macr(I)
         READ(922,*)  WL(I), ACSA_ch3coch3(I), ACSB_ch3coch3(I), ACSC_ch3coch3(I), ACSD_ch3coch3(I)
         READ(923,*)  WL(I), QY_ch3coch3(I)
         READ(924,*)  WL(I), ACS_mek(I), QY_mek(I)
         READ(925,*)  WL(I), ACS_mvk(I)
         READ(926,*)  WL(I), QY1_mvk(I), QY2_mvk(I)
         READ(927,*)  WL(I), ACS_glyox(I)
         READ(928,*)  WL(I), QY1_glyox(I), QY3_glyox(I), QY2_glyox(I)
         READ(929,*)  WL(I), ACS_mglyox(I)
         READ(930,*)  WL(I), QY_mglyox(I)
         READ(931,*)  WL(I), ACS_biacet(I), QY_biacet(I)
         READ(932,*)  WL(I), ACS_ch3ooh(I), QY_ch3ooh(I)
         READ(933,*)  WL(I), ACSA_ch3no3(I), ACSB_ch3no3(I), QY_ch3no3(I)
         READ(934,*)  WL(I), ACSA_c2h5no3(I), ACSB_c2h5no3(I), QY_c2h5no3(I)
         READ(935,*)  WL(I), ACS_n_c3h7no3(I), QY_n_c3h7no3(I)
         READ(936,*)  WL(I), ACSA_i_c3h7no3(I), ACSB_i_c3h7no3(I), QY_i_c3h7no3(I)
         READ(937,*)  WL(I), ACS_t_c4h9no3(I), QY_t_c4h9no3(I)
         READ(938,*)  WL(I), ACS_noa(I), QY1_noa(I), QY2_noa(I)
         READ(939,*)  WL(I), ACS_ho2no2(I), QY1_ho2no2(I), QY2_ho2no2(I)
         READ(940,*)  WL(I), ACSA_n2o5(I), ACSB_n2o5(I), QY1_n2o5(I), QY2_n2o5(I)
         READ(941,*)  WL(I), ACS_test(I)
      ENDDO

      CLOSE(900)
      CLOSE(901)
      CLOSE(902)
      CLOSE(903)
      CLOSE(904)
      CLOSE(905)
      CLOSE(906)
      CLOSE(907)
      CLOSE(908)
      CLOSE(909)
      CLOSE(910)
      CLOSE(911)
      CLOSE(912)
      CLOSE(913)
      CLOSE(914)
      CLOSE(915)
      CLOSE(916)
      CLOSE(917)
      CLOSE(918)
      CLOSE(919)
      CLOSE(920)
      CLOSE(921)
      CLOSE(922)
      CLOSE(923)
      CLOSE(924)
      CLOSE(925)
      CLOSE(926)
      CLOSE(927)
      CLOSE(928)
      CLOSE(929)
      CLOSE(930)
      CLOSE(931)
      CLOSE(932)
      CLOSE(933)
      CLOSE(934)
      CLOSE(935)
      CLOSE(936)
      CLOSE(937)
      CLOSE(938)
      CLOSE(939)
      CLOSE(940)
      CLOSE(941)
    ENDIF ! do only once

    ! Put all J-values to 0 when new time step begins:
    Jpr = 0.0d0

    ! Delta lambda and dimension parameters are given:
    DL  = 5.
    h_planck = 6.62606957D-34 !Js
    c_light = 2.99792458D+17 !nm/s 
    VAA = 1.D-20
    VAB = 1.D-19
    VAC = 1.D-21
    VAD = 1.D-24
    VAAAA = 1.D-20 * 1.D14
    T1 = T - 230.
    AAAA  =  0.332 + 2.5650D-4 * T1 + 1.152D-5 * T1**2 + 2.3130D-8 * T1**3
    BBBB  = -0.575 + 5.5900D-3 * T1 - 1.439D-5 * T1**2 + 3.2700D-8 * T1**3
    CCCC  =  0.466 + 8.8830D-4 * T1 - 3.546D-5 * T1**2 - 3.5190D-7 * T1**3
    LOOO =  308.2 + 4.4871D-2 * T1 + 6.938D-5 * T1**2 - 2.5452D-6 * T1**3

    DO I=1,75
      SAFN_out(I) = SAFN(I) * WL(I)/(h_planck*c_light) / 1.0d4
    ENDDO

    DO I=1,75
      IF (WL(I) .LE. 300) THEN
        QY_o3a = 0.9
      ELSE
        QY_o3a = AAAA * ATAN (BBBB * (REAL(WL(I)) - LOOO)) + CCCC
      ENDIF
      IF (QY_o3a .GT. 0.9) QY_o3a = 0.9
      IF (QY_o3a .LT. 0.)  QY_o3a = 0.
      Jpr(1) = Jpr(1) + ACS_test(I) * QY_o3a * SAFN((WL(I)-280)/5+1) * DL * VAAAA
      IF (Beta .GT. 89.9) THEN
        KULMA = 89.9
        Jpr(2) = 1.23D-3 * EXP(-0.6 / COS(KULMA * 3.14159265 / 180.))
      ELSE
        Jpr(2) = 1.23D-3 * EXP(-0.6 / COS(Beta * 3.14159265 / 180.))
      ENDIF
      IF (Beta .GE. 90.) Jpr(2) = 0.
      
      ! !     ACS_o3 = ACSA_o3(I) * EXP(ACSB_o3(I)/T)
      ! !     PR1 = PR1 + ACS_o3 * QY_o1d(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
      ! !     PR2 = PR2 + ACS_o3 * QY_o3p(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
      ! if (PR3 < -1.0d3 .or. PR3 > 1.0d3 .and. I>70) then
      !   write(*,*) I, PR3, ACS_h2o2(I), QY_h2o2(I), SAFN( (WL(I)-280)/5+1 ), DL, WL(I), (WL(I)-280)/5+1
      ! endif
      ACS_hcho = (ACSA_hcho(I) * VAC) + (ACSB_hcho(I) * VAD * (T-298))
      ACS_ch3coch3 = ACSA_ch3coch3(I) * (1 + (ACSB_ch3coch3(I) * T) + (ACSC_ch3coch3(I) * (T**2)) + (ACSD_ch3coch3(I) * (T**3)))
      ACS_ch3no3 = ACSA_ch3no3(I) * VAA * EXP(ACSB_ch3no3(I) * 1.0d-3 * (T-298))
      ACS_c2h5no3 = ACSA_c2h5no3(I) * EXP(ACSB_c2h5no3(I) * (T-298))
      ACS_i_c3h7no3 = ACSA_i_c3h7no3(I) * EXP(ACSB_i_c3h7no3(I) * (T-298))
      ACS_n2o5 = ACSA_n2o5(I) * 10**(1000*ACSB_n2o5(I)*((1/T)-(1.0d0/298.0d0)))
      
      ! In the old method this is multiplied: SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4 at the end,
      ! but it is just equal to SAFN_out(I).
      ! So Putian modified to a simpler style.
      Jpr(3)  = Jpr(3) + ACS_h2o2(I) * QY_h2o2(I) * SAFN_out(I)
      Jpr(4)  = Jpr(4) + ACS_no2(I) * QY_no2(I) * SAFN_out(I) * VAA
      Jpr(5)  = Jpr(5) + ACS_no3(I) * VAB * QY_no3_no(I) *SAFN_out(I)
      Jpr(6)  = Jpr(6) + ACS_no3(I) * VAB * QY_no3_no2(I) *SAFN_out(I)
      Jpr(7)  = Jpr(7) + ACS_hono(I) * QY_hono(I) *SAFN_out(I)
      Jpr(8)  = Jpr(8) + ACS_hno3(I) * VAA * QY_hno3(I) * SAFN_out(I)
      
      Jpr(11) = Jpr(11) + ACS_hcho * QY_hcho_h_hco(I) * SAFN_out(I)
      Jpr(12) = Jpr(12) + ACS_hcho * QY_hcho_h2_co(I) * SAFN_out(I)
      Jpr(13) = Jpr(13) + ACS_ch3cho(I) * QY_ch3cho(I) * SAFN_out(I)
      Jpr(14) = Jpr(14) + ACS_c2h5cho(I) * QY_c2h5cho(I) * SAFN_out(I)
      Jpr(15) = Jpr(15) + ACS_nc3h7cho(I) * VAC * QY1_nc3h7cho(I) * SAFN_out(I)
      Jpr(16) = Jpr(16) + ACS_nc3h7cho(I) * VAC * QY2_nc3h7cho(I) * SAFN_out(I)
      Jpr(17) = Jpr(17) + ACS_ic3h7cho(I) * VAC * QY_ic3h7cho(I) * SAFN_out(I)
      Jpr(18) = Jpr(18) + ACS_macr(I) * QY1_macr(I) * SAFN_out(I)
      Jpr(19) = Jpr(19) + ACS_macr(I) * QY2_macr(I) * SAFN_out(I)
      Jpr(21) = Jpr(21) + ACS_ch3coch3 * QY_ch3coch3(I) * SAFN_out(I)
      Jpr(22) = Jpr(22) + ACS_mek(I) * QY_mek(I) * SAFN_out(I)
      Jpr(23) = Jpr(23) + ACS_mvk(I) * QY1_mvk(I) * SAFN_out(I)
      Jpr(24) = Jpr(24) + ACS_mvk(I) * QY2_mvk(I) * SAFN_out(I)
      
      Jpr(31) = Jpr(31) + ACS_glyox(I) * VAA * QY1_glyox(I) * SAFN_out(I)
      Jpr(32) = Jpr(32) + ACS_glyox(I) * VAA * QY2_glyox(I) * SAFN_out(I)
      Jpr(33) = Jpr(33) + ACS_glyox(I) * VAA * QY3_glyox(I) * SAFN_out(I)
      Jpr(34) = Jpr(34) + ACS_mglyox(I) * VAA * QY_mglyox(I) * SAFN_out(I)
      Jpr(35) = Jpr(35) + ACS_biacet(I) * VAA * QY_biacet(I) * SAFN_out(I)
      
      Jpr(41) = Jpr(41) + ACS_ch3ooh(I) * VAA * QY_ch3ooh(I) * SAFN_out(I)
      
      Jpr(51) = Jpr(51) + ACS_ch3no3 * QY_ch3no3(I) * SAFN_out(I)
      Jpr(52) = Jpr(52) + ACS_c2h5no3 * QY_c2h5no3(I) * SAFN_out(I)
      Jpr(53) = Jpr(53) + ACS_n_c3h7no3(I) * QY_n_c3h7no3(I) * SAFN_out(I)
      Jpr(54) = Jpr(54) + ACS_i_c3h7no3 * QY_i_c3h7no3(I) * SAFN_out(I)
      Jpr(55) = Jpr(55) + ACS_t_c4h9no3(I) * QY_t_c4h9no3(I) * SAFN_out(I)
      Jpr(56) = Jpr(56) + ACS_noa(I) * QY1_noa(I) * SAFN_out(I)
      Jpr(57) = Jpr(57) + ACS_noa(I) * QY2_noa(I) * SAFN_out(I)
      
      Jpr(67) = Jpr(67) + ACS_ho2no2(I) * QY1_ho2no2(I) * SAFN_out(I) * VAA
      Jpr(68) = Jpr(68) + ACS_ho2no2(I) * QY2_ho2no2(I) * SAFN_out(I) * VAA
      
      Jpr(70) = Jpr(70) + ACS_n2o5(I) * QY1_n2o5(I)  * SAFN_out(I) * VAA
      Jpr(71) = Jpr(71) + ACS_n2o5(I) * QY2_n2o5(I)  * SAFN_out(I) * VAA
    ENDDO

  END SUBROUTINE PHOTOLYSIS_KPP

END MODULE KPP_ROOT_Main
