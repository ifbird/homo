!ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
!
!   CHEMISTRY MODULE
!
!   This module includes:
!       - Conversion of solar irradiance to actinic flux 
!       - Calculation of J-values
!       - Calculation of complex reaction rates 
!       - Calls for KPP that calculates concentrations of the included compounds
!
!   This version is for MCM version 3.2 (J-values & complex reaction rate)
!
!   See version control for comments on changes
!      
!
!oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



MODULE chemistry_mod
  ! Include the following second_*.f90 files that are specific for different chemistries:
  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize

  IMPLICIT NONE

  PRIVATE

  !Public subroutines:
  PUBLIC :: chemistry, kpp_setup

  !Public variables:
  PUBLIC :: NPHOT, NKVALUES
       
  ! REAL(kind=dp) :: T, DVAL(NSPEC)
  REAL(kind=dp) :: RSTATE(20)
  INTEGER :: ICNTRL(20) ! Added by Sampo Smolander
  INTEGER :: i

  !Global variables:
  INTEGER, PARAMETER :: NKVALUES = 46   ! Number of rate coefficients used in KPP. Hand-copied from second_Main.f90
                                        ! NKVALUES = 46 if MCMv3.3 or MCMv3.3.1, NKVALUES = 42 if MCMv3.2

  ! File path for the in- and output values:
  ! CHARACTER(LEN=*), PARAMETER, PRIVATE :: &
  !      filename1 = '../sosa_in', &
  !      filename2 = '../sosa_out'

CONTAINS


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !                                                             
  !   Initialization 
  !                                                             
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE kpp_setup()
    STEPMIN = 0.0d0
    STEPMAX = 0.0d0
    ICNTRL = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /)
    DO i=1,NVAR
       RTOL(i) = 1.0d-4
       ATOL(i) = 1.0d-3
    END DO
    CALL Initialize()
  END SUBROUTINE kpp_setup


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !                                                             
  !   KPP integration
  !                                                             
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE kpp_proceed(C_inout, TSTART_in, TEND_in, TEMP_in, RH_in, PRESS_in, O2_in, N2_in, M_in, H2O_in, RES1_in, RES2_in, J_in, KVALUES_in)

    IMPLICIT NONE

    REAL(kind=dp), INTENT(INOUT) :: C_inout(NSPEC)
    REAL(kind=dp), INTENT(IN   ) :: TSTART_in, TEND_in
    REAL(kind=dp), INTENT(IN   ) :: TEM_in, RH_in, PRESS_in, O2_in, N2_in, M_in, H2O_in, RES1_in, RES2_in
    REAL(kind=dp), INTENT(IN   ) :: J_in(NPHOT),KVALUES_in(NKVALUES)

    REAL(kind=dp) :: T1,T2

    T1 = TSTART_in
    T2 = TEND_in

    TEMP  = TEMP_in
    RH    = RH_in
    PRESS = PRESS_in
    O2    = O2_in
    N2    = N2_in
    M     = M_in
    H2O   = H2O_in
    RES1  = RES1_in
    RES2  = RES2_in
    J     = J_in

    KMT01     = KVALUES_in(1)
    KMT02     = KVALUES_in(2)
    KMT03     = KVALUES_in(3)
    KMT04     = KVALUES_in(4)
    KMT05     = KVALUES_in(5)
    KMT06     = KVALUES_in(6)
    KMT07     = KVALUES_in(7)
    KMT08     = KVALUES_in(8)
    KMT09     = KVALUES_in(9)
    KMT10     = KVALUES_in(10)
    KMT11     = KVALUES_in(11)
    KMT12     = KVALUES_in(12)
    KMT13     = KVALUES_in(13)
    KMT14     = KVALUES_in(14)
    KMT15     = KVALUES_in(15)
    KMT16     = KVALUES_in(16)
    KMT17     = KVALUES_in(17)
    KMT18     = KVALUES_in(18)
    KFPAN     = KVALUES_in(19)
    KBPAN     = KVALUES_in(20)
    KRO2NO    = KVALUES_in(21)
    KRO2HO2   = KVALUES_in(22)
    KAPHO2    = KVALUES_in(23)
    KAPNO     = KVALUES_in(24)
    KRO2NO3   = KVALUES_in(25)
    KNO3AL    = KVALUES_in(26)
    KDEC      = KVALUES_in(27)
    KROPRIM   = KVALUES_in(28)
    KROSEC    = KVALUES_in(29)
    KCH3O2    = KVALUES_in(30)
    K298CH3O2 = KVALUES_in(31)
    K_h2so4   = KVALUES_in(32)
    K_CI      = KVALUES_in(33)
    KMT34     = KVALUES_in(34)
    KMT35     = KVALUES_in(35)
    KMT36     = KVALUES_in(36)
    KMT37     = KVALUES_in(37)
    KMT38     = KVALUES_in(38)
    KMT39     = KVALUES_in(39)
    KMT40     = KVALUES_in(40)
    KMT41     = KVALUES_in(41)
    KMT42     = KVALUES_in(42)
    K16ISOM1  = KVALUES_in(43)
    K15ISOM1  = KVALUES_in(44)
    K14ISOM1  = KVALUES_in(45)
    KBPPN     = KVALUES_in(46)

    ! Update_RCONST and INTEGRATE will operate on C
    C = C_inout
    CALL Update_RCONST()
    CALL INTEGRATE( TIN = T1, TOUT = T2, RSTATUS_U = RSTATE, ICNTRL_U = ICNTRL )
    C_inout = C
  END SUBROUTINE kpp_proceed


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !                                                             
  !   CHEMIRTY                                               
  !                                                             
  !   Calculates the chemistry for each layer by using KPP                 
  !                                                             
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE chemistry(CONS, TIME_kpp, END_kpp, tmcontr, year_int, month_int, nxodrad, Tp, Pr, ZSD_in, SUN_PAR_LAYER, &
       RES1, RES2, q_H2O, H2O, RO2_out, J_values, K_values, chem_glo, STATION, Aground, I_ACF_out)
    

    IMPLICIT NONE

    INTEGER                                         :: i20, j20, J, I, j1, loc
    INTEGER, SAVE                                   :: FIRST_CALL = 1
    INTEGER, INTENT(in)                             :: year_int,           & ! What year your simulation is done for
                                                       month_int,          & ! What month your simulation is done for
                                                       nxodrad               ! Time number for every half hour

    CHARACTER(len=3)                                :: STATION,            & ! What station this is
                                                       month                 ! What month this is
    CHARACTER(len=4)                                :: year                  ! what year this is
    CHARACTER(len=3), SAVE                          :: month_old = 'M00'     ! So that in the very beginning of the program, month_old initial 
                                                                             ! value matches no real month (real month values are 'M01' ... 'M12')

    REAL(kind=dp), PARAMETER                        :: Avog   = 6.0223d23, & ! Avogrado's constant in molecules / mol
                                                       RG     = 8.314d0,     & ! Gas constant in J / mol / K
                                                       M_Air  = 28.97d0,     & ! Average molar mass of air in g / mol
                                                       M_H2O  = 18.0153d0      ! Molar mass of water in g / mol

    REAL(kind=dp)                                   :: ZSD,                & ! Solar zenith angle
                                                       Air,                & ! Concentration of air molecules
                                                       q_H2O,              &
                                                       RH,                 & ! Relative humidity 
                                                       O2,                 & ! Concentration of O2 in molecule / cm3
                                                       N2,                 & ! Concentration of N2 in molecule / cm3
                                                       MO2N2,              & ! Sum of O2 and N2 concentration molecule / cm3
                                                       H2,                 & ! Concentration of H2 in molecule / cm3
                                                       H2O,                & ! Concentration of H2O in molecule / cm3
                                                       TEMP_KPP,           & ! Temperature used in KPP. Unit: K
                                                       Pr_KPP,             & ! Pressure used in KPP
                                                       ZSD_n                 ! Parameter to change ZSD from REAL to INTEGER

    REAL(kind=dp), INTENT(inout)                    :: CONS(NSPEC)           ! Gas concentrations

    REAL(kind=dp), INTENT(in)                       :: TIME_kpp,           & ! Time step in KPP
                                                       END_kpp,            & ! End time in KPP
                                                       RES1,               & ! 'reaction rate' for how fast H2SO4 goes to particle phase
                                                       RES2,               & ! 'reaction rate' for how fast HNO3 goes to particle phase
                                                       tmcontr,            & ! Time in seconds
                                                       Tp,                 & ! temperature at this level of atmosphere
                                                       Pr,                 & ! pressure at this level of atmosphere
                                                       ZSD_in,             & ! Solar zenith angle
                                                       SUN_PAR_LAYER,      & ! Penetration of radiation trough the canopy. Scales the radiation, so there is less radiation in the floor of the canopy
                                                       chem_glo,           & ! incoming global radiation (W/m2)
                                                       Aground               ! ground albedo from measurements (not sure for what wl)
  
    REAL(kind=dp), INTENT(out)                      :: RO2_out               ! Concentration of peroxy radical    

    REAL(kind=dp), INTENT(out), DIMENSION(NPHOT)    :: J_values              ! J_values, so 'rate constants' for photolysis reactions
    REAL(kind=dp), INTENT(out), DIMENSION(NKVALUES) :: K_values              ! K_values; complex rate constants
    REAL(kind=dp), INTENT(out), DIMENSION(75)       :: I_ACF_out             ! Calculated actinic flux. Unit: photon * cm^-2 * s^-1 * nm^-1

    REAL(kind=dp), DIMENSION(47,1488), SAVE         :: swr_hyy               !Measured spectral irradiance. Unit: mW/(m^2*nm)
    REAL(kind=dp), DIMENSION(84), SAVE              :: swr_distribution
    REAL(kind=dp), DIMENSION(75)                    :: I_ACF                 ! Calculated actinic flux. Different unit than I_ACF_out
    REAL(kind=dp), DIMENSION(45,90), SAVE           :: E_o_E                 ! wl and ZSD dependent conversion factor from irradiance to actinic flux

!----------------------------------------------------------------------------------------------------------------------

    WRITE(year,'(I4)') year_int
    IF (month_int < 10) WRITE(month,'("M0",I1)') month_int
    IF (month_int >= 10) WRITE(month,'("M",I2)') month_int

    IF (STATION .EQ. 'MAN') THEN
       IF (FIRST_CALL .EQ. 1) THEN
          FIRST_CALL = 0
          OPEN(unit=4,file = ''//filename1//'/General/Manitou/swr_distribution.txt',status='old')
          READ(4,*) (swr_distribution(i20),i20=1,84)
          CLOSE(unit=4)
       ENDIF
    ENDIF

    !IF (STATION .EQ. 'HYY') THEN  !Currently we are not using this because of
    !non-calibrated spectral radiation data
    !   IF (month_old /= month) THEN ! so that this is run once only when the month has changed
    !      month_old = month
    !      ! Input of measured spectral radiation data (every half-hour data)
    !      OPEN(unit=4,file = ''//filename1//'/'//STATION//'/Year/'//YEAR//'/Data_input/'//MONTH//'/hyy_swr.txt',status='old')
    !      READ(4,*)((swr_hyy(i20,j20),i20=1,47),j20=1,1488)
    !      CLOSE(unit=4)
    !
    !      ! wl and ZSD dependent conversion factor from irradiance to actinic flux
    !      OPEN(unit=4,file=''//filename1//'/General/Hyytiala/solar_factor.txt',status='old')
    !      READ(4,*)((E_o_E(i20,j20),i20=1,45),j20=1,90)
    !      CLOSE(unit=4)
    !   ENDIF
    !ENDIF
    !***** debug for E_o_E *****!
    ! write(*,*) E_o_E(1,:)
    ! do i=1,45
    !   write(*,*) i, SUM(E_o_E(i,:))
    ! end do
    ! write(*,*) SUM(E_o_E)
    !***** end debug *****!

    IF (STATION .EQ. 'HYY') THEN 
       IF (FIRST_CALL .EQ. 1) THEN
          FIRST_CALL = 0
          OPEN(unit=4,file=''//filename1//'/General/swr_distribution.txt',status='old')
          READ(4,*) (swr_distribution(i20),i20=1,84)
          CLOSE(unit=4)
       ENDIF
    ENDIF

    ! Calculation of Air, H2O, O2, N2, H2 and MO2N2 in molecule / cm3
    ZSD      = ZSD_in
    Air      = Pr / Tp / RG * Avog / 1.0d4 
    O2       = 0.209460d0 * Air
    N2       = 0.780840d0 * Air
    H2       = 0.5d0 * Air / 1.0d6
    MO2N2    = O2 + N2
    H2O      = q_H2O / M_H2O * Avog / 1.0d3
    !Aground  = 0.1 ! should normally come from measurements if available
    TEMP_KPP = Tp    ! because TEMP_KPP has type KIND=dp
    Pr_KPP   = Pr    ! same as for temp
    RH       = 50   ! not in use at this time


    CALL ratecons(Tp, O2, N2, MO2N2, H2O, K_values, TIME_kpp) ! Calculates the rates for complex gas-reactions 


    ! Calculated actinic flux from spectral data measured in Hyytiala 
    IF (zsd .eq. 0.) THEN
      zsd = 90.
    ENDIF
 
    if (station .eq. 'MAN') then
      do j = 1, 75
        I_ACF(j) = chem_glo * swr_distribution(j) * (1 + 2 * Aground * COS(ZSD*3.14/180) + Aground)
      enddo
    elseif (station .eq. 'HYY') then
      !ZSD_n = ANINT(ZSD)
      !LOC = INT(ZSD_n)
      !DO J = 1,45 ! For solar radiation: measured data is read in for 280-500 nm.
      !  I_ACF(J) = E_o_E(J,LOC) * SUN_PAR_LAYER / 1000.0d0 / 5.0d0 * &
      !             (swr_hyy(J+1,nxodrad) + (tmcontr-(nxodrad-1)*1800.0d0) * (swr_hyy(J+1,nxodrad+1) - swr_hyy(J+1,nxodrad))/1800.0d0)
      !ENDDO
      !DO J = 46,75 ! Solar radiation for 500 nm is used for 505 - 650 nm. 
      !  I_ACF(J) = I_ACF(45)
      !ENDDO       
       do j = 1, 75
          I_ACF(j) = chem_glo * swr_distribution(j) * (1 + 2 * Aground*COS(ZSD*3.14/180) + Aground)
       enddo
    endif
    

    CALL PHOTOLYSIS_KPP   (tmcontr, ZSD, I_ACF, I_ACF_out, TP,                                                           &
         J_values(1), J_values(2), J_values(3), J_values(4), J_values(5), J_values(6), J_values(7),                      &
         J_values(8), J_values(11), J_values(12), J_values(13), J_values(14), J_values(15), J_values(16), J_values(17),  &
         J_values(18), J_values(19), J_values(21), J_values(22), J_values(23), J_values(24), J_values(31), J_values(32), &
         J_values(33), J_values(34), J_values(35), J_values(41), J_values(51), J_values(52), J_values(53), J_values(54), &
         J_values(55), J_values(56), J_values(57), J_values(67), J_values(68), J_values(70), J_values(71)) ! Calculates J-values



    IF (ZSD .GE. 90) THEN
       DO J = 1, NPHOT
          J_values(J) = 0.
       ENDDO
    ENDIF

    ! Some safety check
    DO J1 = 1,NPHOT
       IF (J_values(J1) < -1e9 .OR. J_values(J1) > 1e9 ) then
          WRITE(*,*) 'Note by Sampo: J_value(J1) has bad value.'
          WRITE(*,*) 'J1, J_value(J1) = ', J1, J_values(J1)
          STOP
       ENDIF
    ENDDO


    CALL KPP_Proceed(CONS, TIME_kpp, END_kpp, TEMP_KPP, RH, Pr_KPP, O2, N2, MO2N2, H2O, RES1, RES2, J_values, K_values)

    RO2_out = RO2

  END SUBROUTINE chemistry


  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !                                                                           
  !   Subroutine for the photodissociation reactions  
  !                                                                         
  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo


  SUBROUTINE PHOTOLYSIS_KPP(tmcontr, Beta, SAFN, SAFN_out, T,       &
       PR1, PR2, PR3, PR4, PR5, PR6, PR7, PR8, PR11,                &
       PR12, PR13, PR14, PR15, PR16, PR17, PR18, PR19, PR21, PR22,  &
       PR23, PR24, PR31, PR32, PR33, PR34, PR35, PR41, PR51, PR52,  &
       PR53, PR54, PR55, PR56, PR57, PR67, PR68, PR70, PR71)

    IMPLICIT NONE

    INTEGER ::                      I

    INTEGER, DIMENSION(75), SAVE :: WL !wavelength

    REAL(kind=dp), DIMENSION(75) :: SAFN, SAFN_out ! Actinic flux in two different units

    REAL(kind=dp) :: T,                                                                 & ! temperature in K
         tmcontr,                                                                       & ! time in s
         Beta,                                                                          & ! solar zenit angle
         DL,                                                                            & ! Delta lambda (so wl band)
         VAA, VAB, VAC, VAD, h_planck, c_light,                                         & ! Dimension parameter
         PR1, PR2, PR3, PR4, PR5, PR6, PR7, PR8, PR11, PR12, PR13, PR14, PR15,          & ! PRX : photolysis rate 
         PR16, PR17, PR18, PR19, PR21, PR22, PR23, PR24, PR31, PR32, PR33, PR34,        &
         PR35, PR41, PR51, PR52, PR53, PR54, PR55, PR56, PR57, PR67, PR68, PR70, PR71,  &
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
    PR1 = 0.
    PR2 = 0.
    PR3 = 0.
    PR4 = 0.
    PR5 = 0.
    PR6 = 0.
    PR7 = 0.
    PR8 = 0.
    PR11 = 0.
    PR12 = 0.
    PR13 = 0.
    PR14 = 0.
    PR15 = 0.
    PR16 = 0.
    PR17 = 0.
    PR18 = 0.
    PR19 = 0.
    PR21 = 0.
    PR22 = 0.
    PR23 = 0.
    PR24 = 0.
    PR31 = 0.
    PR32 = 0.
    PR33 = 0.
    PR34 = 0.
    PR35 = 0.
    PR41 = 0.
    PR51 = 0.
    PR52 = 0.
    PR53 = 0.
    PR54 = 0.
    PR55 = 0.
    PR56 = 0.
    PR57 = 0.
    PR67 = 0.
    PR68 = 0.
    PR70 = 0.
    PR71 = 0.

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
      SAFN_out(I) = SAFN(I) * WL(I)/(h_planck*c_light) * 1/10**4 
   ENDDO

   DO I=1,75
         IF (WL(I) .LE. 300) THEN
            QY_o3a = 0.9
         ELSE
            QY_o3a = AAAA * ATAN (BBBB * (REAL(WL(I)) - LOOO)) + CCCC
         ENDIF
         IF (QY_o3a .GT. 0.9) QY_o3a = 0.9
         IF (QY_o3a .LT. 0.)  QY_o3a = 0.
         PR1 = PR1 + ACS_test(I) * QY_o3a * SAFN((WL(I)-280)/5+1) * DL * VAAAA
         IF (Beta .GT. 89.9) THEN
         KULMA = 89.9
         PR2 = 1.23D-3 * EXP(-0.6 / COS(KULMA * 3.14159265 / 180.))
         ELSE
         PR2 = 1.23D-3 * EXP(-0.6 / COS(Beta * 3.14159265 / 180.))
         ENDIF
         IF (Beta .GE. 90.) PR2 = 0.

   ! !     ACS_o3 = ACSA_o3(I) * EXP(ACSB_o3(I)/T)
   ! !     PR1 = PR1 + ACS_o3 * QY_o1d(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
   ! !     PR2 = PR2 + ACS_o3 * QY_o3p(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
      ! if (PR3 < -1.0d3 .or. PR3 > 1.0d3 .and. I>70) then
      !   write(*,*) I, PR3, ACS_h2o2(I), QY_h2o2(I), SAFN( (WL(I)-280)/5+1 ), DL, WL(I), (WL(I)-280)/5+1
      ! endif
         PR3 = PR3 + ACS_h2o2(I) * QY_h2o2(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR4 = PR4 + ACS_no2(I) * QY_no2(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4 * VAA
         PR5 = PR5 + ACS_no3(I) * VAB * QY_no3_no(I) *SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR6 = PR6 + ACS_no3(I) * VAB * QY_no3_no2(I) *SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR7 = PR7 + ACS_hono(I) * QY_hono(I) *SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR8 = PR8 + ACS_hno3(I) * VAA * QY_hno3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         ACS_hcho = (ACSA_hcho(I) * VAC) + (ACSB_hcho(I) * VAD * (T-298))
         PR11 = PR11 + ACS_hcho * QY_hcho_h_hco(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR12 = PR12 + ACS_hcho * QY_hcho_h2_co(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR13 = PR13 + ACS_ch3cho(I) * QY_ch3cho(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR14 = PR14 + ACS_c2h5cho(I) * QY_c2h5cho(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR15 = PR15 + ACS_nc3h7cho(I) * VAC * QY1_nc3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR16 = PR16 + ACS_nc3h7cho(I) * VAC * QY2_nc3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR17 = PR17 + ACS_ic3h7cho(I) * VAC * QY_ic3h7cho(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR18 = PR18 + ACS_macr(I) * QY1_macr(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR19 = PR19 + ACS_macr(I) * QY2_macr(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         ACS_ch3coch3 = ACSA_ch3coch3(I) * (1 + (ACSB_ch3coch3(I) * T) + (ACSC_ch3coch3(I) * (T**2)) + (ACSD_ch3coch3(I) * (T**3)))
         PR21 = PR21 + ACS_ch3coch3 * QY_ch3coch3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR22 = PR22 + ACS_mek(I) * QY_mek(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR23 = PR23 + ACS_mvk(I) * QY1_mvk(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR24 = PR24 + ACS_mvk(I) * QY2_mvk(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR31 = PR31 + ACS_glyox(I) * VAA * QY1_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR32 = PR32 + ACS_glyox(I) * VAA * QY2_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR33 = PR33 + ACS_glyox(I) * VAA * QY3_glyox(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR34 = PR34 + ACS_mglyox(I) * VAA * QY_mglyox(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR35 = PR35 + ACS_biacet(I) * VAA * QY_biacet(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR41 = PR41 + ACS_ch3ooh(I) * VAA * QY_ch3ooh(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         ACS_ch3no3 = ACSA_ch3no3(I) * VAA * EXP(ACSB_ch3no3(I) * 1.E-3 * (T-298))
         PR51 = PR51 + ACS_ch3no3 * QY_ch3no3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         ACS_c2h5no3 = ACSA_c2h5no3(I) * EXP(ACSB_c2h5no3(I) * (T-298))
         PR52 = PR52 + ACS_c2h5no3 * QY_c2h5no3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR53 = PR53 + ACS_n_c3h7no3(I) * QY_n_c3h7no3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         ACS_i_c3h7no3 = ACSA_i_c3h7no3(I) * EXP(ACSB_i_c3h7no3(I) * (T-298))
         PR54 = PR54 + ACS_i_c3h7no3 * QY_i_c3h7no3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR55 = PR55 + ACS_t_c4h9no3(I) * QY_t_c4h9no3(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR56 = PR56 + ACS_noa(I) * QY1_noa(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR57 = PR57 + ACS_noa(I) * QY2_noa(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4
         PR67 = PR67 + ACS_ho2no2(I) * QY1_ho2no2(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4 * VAA
         PR68 = PR68 + ACS_ho2no2(I) * QY2_ho2no2(I) * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4 * VAA
         ACS_n2o5 = ACSA_n2o5(I) * 10**(1000*ACSB_n2o5(I)*((1/T)-(1/298)))
         PR70 = PR70 + ACS_n2o5(I) * QY1_n2o5(I)  * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4 * VAA
         PR71 = PR71 + ACS_n2o5(I) * QY2_n2o5(I)  * SAFN((WL(I)-280)/5+1) * DL * WL(I)/(h_planck*c_light) * 1/10**4 * VAA
    ENDDO

  END SUBROUTINE PHOTOLYSIS_KPP



  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !                                                                         
  !   ratecons                                     
  !                                                                         
  !   Calculation of simple and complex rate coefficients 
  !                                                                         
  !oooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo



  SUBROUTINE ratecons(Temp, O2, N2, M, H2O, RCON)
    REAL(kind=dp), INTENT(IN) :: Temp, O2, N2, M, H2O
    REAL(kind=dp), DIMENSION(NKVALUES), INTENT(OUT) :: RCON

    ! Intermediate variables
    REAL(kind=dp) :: K0, KI, KR, FC, K1, K2, K3, K4, NC, FCC, avo, m_h2so4, rhov, amu, Kb, pil
    REAL(kind=dp) :: T3
    ! REAL(kind=dp) :: KRO2NO, KRO2HO2, KAPHO2, KAPNO, KRO2NO3, KNO3AL, KDEC, KROPRIM, KROSEC, KCH3O2, K298CH3O2, K_h2so4

    ! Initiate
    RCON = 0.0d0
    T3 = Temp/300

    !====================================================================!
    ! Complex rate coefficients
    !
    ! Reference: http://mcm.leeds.ac.uk/MCMv3.3.1/parameters/complex.htt
    !====================================================================!

    ! O(1P) + NO -> NO2,   KMT01
    K0        =  1.0d-31 * T3**(-1.6d0) * M
    KI        =  5.0d-11 * T3**(-0.3d0)
    FC        =  0.85d0
    RCON(1)   = f_RCON(K0, KI, FC)

    ! O(1P) + NO2 -> NO3,  KMT02
    K0        =  1.3d-31 * T3**(-1.5d0) * M
    KI        =  2.3d-11 * T3**0.24d0
    FC        =  0.6d0
    RCON(2)   = f_RCON(K0, KI, FC)

    ! NO2 + NO3 -> N2O5,   KMT03
    K0        =  3.6d-30 * T3**(-4.1d0) * M
    KI        =  1.9d-12 * T3**0.2d0
    FC        =  0.35d0
    RCON(3)   = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(3)   = (K0*KI) * FCC / (K0+KI)

    ! N2O5 -> NO2 + NO3,    KMT04
    K0        =  1.3d-3  * T3**(-3.5d0) * EXP(-11000/Temp) * M
    KI        =  9.7d+14 * T3**0.1d0 * EXP(-11080/Temp)
    FC        =  0.35d0
    RCON(4)   = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(4)   = (K0*KI) * FCC / (K0+KI)

    ! OH + CO -> HO2 + CO2,  KMT05
    RCON(5)   =  1.44D-13*(1+(M/4.2D+19))

    !HO2 + HO2        ->  H2O2 + O2,   KMT06
    RCON(6)   =  1 + (1.4d-21 * EXP(2200/Temp) * H2O)

    ! OH + NO -> HONO,    KMT07
    K0        =  7.40d-31 * T3**(-2.4d0) * M
    KI        =  3.30d-11 * T3**(-0.3d0)
    FC        =  0.81d0
    RCON(7)   = f_RCON(K0, KI, FC)
    ! FC        =  EXP(-Temp / 1420.)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(7)   = (K0*KI) * FCC / (K0+KI)

    ! OH + NO2 ->  HNO3,    KMT08
    K0        =  3.2d-30 * T3**(-4.5d0) * M
    KI        =  3.0d-11
    FC        =  0.41d0
    RCON(8)   = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(8)   = (K0*KI) * FCC / (K0+KI)

    ! HO2 + NO2 -> HO2NO2,    KMT09 
    K0        =  1.4d-31 * T3**(-3.1d0) * M
    KI        =  4.0d-12
    FC        =  0.4d0
    RCON(9)   = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(9)   = (K0*KI) * FCC / (K0+KI)

    ! HO2NO2 -> HO2 + NO2,    KMT10
    K0        =  4.1d-5 * EXP(-10650 / Temp) * M
    KI        =  6.0d+15 * EXP(-11170 / Temp)
    FC        = 0.4d0
    RCON(10)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(10)   = (K0*KI) * FCC / (K0+KI)

    ! OH + HNO3 -> NO3 + H2O,    KMT11
    K1        =  2.4d-14 * EXP(460/Temp)
    K3        =  6.5d-34 * EXP(1335/Temp)
    K4        =  2.7d-17 * EXP(2199/Temp)
    K2        =  (K3 * M) / (1 + (K3*M/K4))
    RCON(11)  =  K1 + K2

    ! OH + SO2 -> HSO3,    KMT12
    K0        =  2.5d-31 * T3**(-2.6d0) * M
    KI        =  2.0d-12
    FC        =  0.53d0
    RCON(12)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(12)   = (K0*KI* FCC) / (K0+KI)

    !Some old way

    !K0        =  4.0E-31 * (Temp/300.)**(-3.3) * M
    !KI        =  2.0E-12
    !FC        =  0.45
    !RCON(12)  =  (K0 / (1. + K0/KI)) * FC**(1 / (1 + (LOG10(K0 / KI))**2))

    !S&P / Aitkinson
    !K0        =  3.0E-31 * (Temp/300.)**(-3.3) * M
    !KI        =  1.5E-12
    !FC        =  0.6
    !RCON(12)  =  (K0 / (1. + K0/KI)) * FC**(1 / (1 + (LOG10(K0 / KI))**2))

    ! CH3O2 + NO2 -> CH3O2NO2,   KMT13
    K0        =  2.50d-30 * T3**(-5.5d0) * M
    KI        =  1.8d-11
    FC        =  0.36d0
    RCON(13)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(13)   = (K0*KI) * FCC / (K0+KI)

    ! CH3O2NO2 -> CH3O2 + NO2,    KMT14
    K0        =  9.00d-5 * EXP(-9690/Temp) * M
    KI        =  1.10d+16 * EXP(-10560/Temp)
    FC        =  0.36d0
    RCON(14)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(14)   = (K0*KI) * FCC / (K0+KI)

    ! OH + C2H4 -> C2H4OH,    KMT15
    K0        =  8.6d-29 * T3**(-3.1d0) * M
    KI        =  9.00d-12* T3**(-0.85d0)
    FC        =  0.48d0
    RCON(15)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(15)   = (K0*KI) * FCC / (K0+KI)

    ! OH + C3H6 -> C3H6OH,    KMT16
    K0        =  8.00d-27 * T3**(-3.5d0) * M
    KI        =  3.00d-11 * T3**(-1)
    FC        =  0.5d0
    RCON(16)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(16)   = (K0*KI) * FCC / (K0+KI)

    ! OH + C2H2 -> C2H2OH,    KMT17
    K0        =  5.0d-30 * T3**(-1.5d0) * M
    KI        =  1.0d-12
    FC        =  0.17d0*EXP(-51/Temp) + EXP(-Temp/204)
    RCON(17)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(17)   = (K0*KI*FCC) / (K0+KI)

    ! OH + CH3SCH3 -> 4??   KMT18
    RCON(18)  =  9.5d-39 * O2 * EXP(5270/Temp) / (1 + 7.5d-29 * O2 * EXP(5610/Temp))

    ! KFPAN - MCM v3.2
    !K0        =  2.70E-28 * EXP(Temp/300)**(-7.1) * M
    !KI        =  1.20E-11 * EXP(Temp/300)**(-0.9)
    !FC        =  0.3
    !NC        =  0.75 - 1.27 * (LOG10(FC))
    !FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    !RCON(19)  =  (K0 * KI) * FCC/(K0+KI)

    ! CH3C(O)O2 + NO2 + M -> CH3C(O)O2NO2 + M,  KFPAN - MCM v3.3
    K0        =  3.28d-28 * M* T3**(-6.87d0)
    KI        =  1.125d-11 * T3**(-1.105d0)
    FC        =  0.3d0
    RCON(19)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(19)  =  (K0 * KI) * FCC/(K0+KI)

    ! KBPAN - MCM v3.2
    !K0        =  4.90E-03 * EXP(-12100/Temp) * M
    !KI        =  5.40E+16 * EXP(-13830/Temp)   !Previously there was a mistake so that it was 5.4E-D16
    !FC        =  0.3
    !NC        =  0.75 - 1.27 * (LOG10(FC))
    !FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    !RCON(20)  =  (K0*KI) * FCC / (K0+KI)

    ! CH3C(O)O2NO2 + M -> CH3C(O)O2 + NO2 + M, KBPAN - MCM v3.3
    K0        = 1.10d-05 * M* EXP(-10100/Temp)
    KI        = 1.90d17 * EXP(-14100/Temp)
    FC        =  0.3d0
    RCON(20)  = f_RCON(K0, KI, FC)
    ! NC        =  0.75 - 1.27 * (LOG10(FC))
    ! FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(20)  =  (K0*KI) * FCC / (K0+KI)


    !Here comes the simple rate coefficients:
    ! KRO2NO
    RCON(21) = 2.7D-12*EXP(360/Temp)

    ! KRO2HO2
    RCON(22) = 2.91D-13*EXP(1300/Temp)

    ! KAPHO2
    RCON(23) = 5.2D-13*EXP(980/Temp)

    ! KAPNO
    RCON(24) = 7.5D-12*EXP(290/Temp)

    ! KRO2NO3
    RCON(25) = 2.3D-12

    ! KNO3AL
    RCON(26) = 1.4D-12*EXP(-1860/Temp)

    ! KDEC
    RCON(27) = 1.00D+06

    ! KROPRIM
    RCON(28) = 2.50D-14*EXP(-300/Temp)

    ! KROSEC
    RCON(29) = 2.50D-14*EXP(-300/Temp)

    ! KCH3O2
    RCON(30) = 1.03D-13*EXP(365/Temp)

    ! K298CH3O2
    RCON(31) = 3.5D-13

    !K_h2so4, collision rate between two H2SO4 molecules
    avo  = 6.0221D23 ! avogadro no.,   /mol
    m_h2so4    = 98.0775 ! molar mass,  g/mol
    rhov = 1830. ! vapor molecule density,  kg/m^3
    amu  = 1.6605D-27 ! kg
    Kb   = 1.380658D-23 ! J/K
    pil = 3.14
    FC = amu * (m_h2so4 / 2.)
    FCC = ((m_h2so4 / avo) * (1D-03 / rhov) * (6. / pil))**(1./3.)
    NC = SQRT(8. * Kb * Temp / (pil * FC))
    RCON(32) = pil * FCC**2. * NC * 10.**6.


    !Criegee intermediate reaction rate: K_CI
    !RCON(33) = 7.00D-14

!!!Now come the rate coefficients from Aitkinson 2004
    ! 2OH-> H2O2,   KMT34
    K0        =  6.9D-31*(Temp/300.)**(-0.8)*M
    KI        =  2.6D-11
    FC        =  0.50
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(34)   = (K0*KI) * FCC / (K0+KI)

    ! NO + NO2 -> N2O3,   KMT35
    K0        =  3.1D-34*(Temp/300.)**(-7.7)*M
    KI        =  7.9D-12*(Temp/300)**1.4
    FC        =  0.6
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(35)   = (K0*KI) * FCC / (K0+KI)

    ! N2O3 + M -> NO + NO2 + M,   KMT36
    K0        =  1.9D-7*(Temp/300.)**(-8.7)*EXP(-4880/Temp)*M
    KI        =  4.7D15*(Temp/300)**(0.4)*EXP(-4880/Temp)
    FC        =  0.6
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(36)   = (K0*KI) * FCC / (K0+KI)

    ! 2NO2 + M -> N2O4 + M,   KMT37
    K0        =  1.4D-33*(Temp/300.)**(-3.8)*M
    KI        =  1.D-12
    FC        =  0.4
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(37)   = (K0*KI) * FCC / (K0+KI)

    ! N2O4 + M -> 2NO2 + M,   KMT38
    K0        =  1.3D-5*(Temp/300.)**(-3.8)*EXP(-6400/Temp)*M
    KI        =  1.15D16*EXP(-6400/Temp)
    FC        =  0.4
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(38)   = (K0*KI) * FCC / (K0+KI)

    ! OH + CS2 + M -> HOCS2 + M,   KMT39
    K0        =  8.D-31*M
    KI        =  8.D-12
    FC        =  0.8
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(39)   = (K0*KI) * FCC / (K0+KI)

    ! HOCS2 + M -> OH + CS2 + M,   KMT40
    K0        =  1.6D-6*EXP(-5160/Temp)*M
    KI        =  1.6D13*EXP(-5160/Temp)
    FC        =  0.8
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(40)   = (K0*KI) * FCC / (K0+KI)

    ! HS + NO + M -> HSNO + M,   KMT41
    K0        =  2.4D-31*(Temp/300)**(-2.5)*M
    KI        =  2.7D-11
    FC        =  0.6
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(41)   = (K0*KI) * FCC / (K0+KI)

    ! CH3S + NO + M -> CH3SNO + M,   KMT42
    K0        =  3.3D-29*(Temp/300)**(-4)*M
    KI        =  4.D-11
    FC        =  0.54
    NC        =  0.75 - 1.27 * (LOG10(FC))
    FCC       =  10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    RCON(42)   = (K0*KI) * FCC / (K0+KI)

    !More simple rates for MCM v3.3
    !K16ISOM1
    RCON(43) = 4.60D9*EXP(-8380/Temp)*EXP(1.00D8/Temp**3)

    !K15ISOM1
    RCON(44) = 1.50D11*EXP(-9750/Temp)

    !K14ISOM1
    RCON(45) = 3.00D7*EXP(-5300/Temp)

    ! KBPPN
    K0        = 1.7d-03 * EXP(-11280/Temp) * M
    KI        = 8.3d+16 * EXP(-13940/Temp)
    FC        = 0.36d0
    RCON(46)  = f_RCON(K0, KI, FC)
    ! NC        = 0.75 - 1.27 * (LOG10(FC))
    ! FCC       = 10**(LOG10(FC)/(1+(LOG10(K0/KI)/NC)**2))
    ! RCON(46)  = (K0*KI) * FCC / (K0+KI)


    CONTAINS

      FUNCTION f_NC(FC)
        REAL(kind=dp) :: FC
        REAL(kind=dp) :: f_NC

        f_NC = 0.75d0 - 1.27d0 * LOG10(FC)
      END FUNCTION f_NC


      FUNCTION f_FCC(K0, KI, FC)
        REAL(kind=dp) :: K0, KI, FC
        REAL(kind=dp) :: f_FCC

        f_FCC = 10**( LOG10(FC) / (1 + (LOG10(K0/KI)/f_NC(FC))**2) )
      END FUNCTION f_FCC


      FUNCTION f_RCON(K0, KI, FC)
        REAL(kind=dp) :: K0, KI, FC
        REAL(kind=dp) :: f_RCON

        f_RCON = (K0 * KI) * f_FCC(K0, KI, FC) / (K0 + KI)
      END FUNCTION f_RCON

  END SUBROUTINE ratecons

END MODULE Chemistry_Mod
