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
  USE mcm_Main
  USE mcm_Model
  USE mcm_Initialize, ONLY: Initialize

  IMPLICIT NONE

  PRIVATE

  !Public subroutines:
  PUBLIC :: chemistry, kpp_setup

  !Public variables:
  PUBLIC :: NPHOT
       
  ! REAL(kind=dp) :: T, DVAL(NSPEC)
  REAL(kind=dp) :: RSTATE(20)
  INTEGER :: ICNTRL(20) ! Added by Sampo Smolander
  INTEGER :: i

  !Global variables:


CONTAINS


  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo
  !                                                             
  !   CHEMIRTY                                               
  !                                                             
  !   Calculates the chemistry for each layer by using KPP                 
  !                                                             
  !ooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooooo

  SUBROUTINE chemistry(C_inout, TSTART_in, TEND_in, TEMP_in, RH_in, PRESS_in, RGlob_in, zsd_in, Aground_in, &
       O2_in, N2_in, M_in, H2O_in, RES1_in, RES2_in, I_ACF_in, RO2_out, J_values, I_ACF_out)
    

    IMPLICIT NONE

    REAL(kind=dp), INTENT(inout) :: C_inout(NSPEC)  ! [molec cm-3], Gas concentrations

    REAL(kind=dp), INTENT(IN) :: TSTART_in,  &  ! Time step in KPP
                                 TEND_in,    &  ! End time in KPP
                                 TEMP_in,    &  ! [K], temperature at this level of atmosphere
                                 RH_in,      &  ! [-], relative humidity
                                 PRESS_in,   &  ! [Pa], pressure at this level of atmosphere
                                 RGlob_in,   &  ! incoming global radiation (W/m2)
                                 zsd_in,     &  ! [??], Solar zenith angle
                                 Aground_in     ! ground albedo from measurements (not sure for what wl)

    REAL(kind=dp), INTENT(IN) :: O2_in,    &  ! [molec cm-3]
                                 N2_in,    &
                                 M_in,     &
                                 H2O_in,   &
                                 RES1_in,  &  ! 'reaction rate' for how fast H2SO4 goes to particle phase
                                 RES2_in      ! 'reaction rate' for how fast HNO3 goes to particle phase
    REAL(kind=dp), INTENT(IN), DIMENSION(75) :: I_ACF_in  ! Calculated actinic flux. Different unit than I_ACF_out

    REAL(kind=dp), INTENT(OUT)                      :: RO2_out    ! [molec cm-3], Concentration of peroxy radical    
    REAL(kind=dp), INTENT(OUT), DIMENSION(NPHOT)    :: J_values   ! J_values, so 'rate constants' for photolysis reactions
    REAL(kind=dp), INTENT(OUT), DIMENSION(75)       :: I_ACF_out  ! Calculated actinic flux. Unit: photon * cm^-2 * s^-1 * nm^-1

    REAL(kind=dp), PARAMETER  :: Avog   = 6.0223d23, & ! Avogrado's constant in molecules / mol
                                 RG     = 8.314d0,     & ! Gas constant in J / mol / K
                                 M_Air  = 28.97d0,     & ! Average molar mass of air in g / mol
                                 M_H2O  = 18.0153d0      ! Molar mass of water in g / mol

    REAL(kind=dp)             :: zsd


    ! Calculated actinic flux from spectral data measured in Hyytiala 
    zsd = zsd_in
    IF (zsd == 0.0d0) THEN
      zsd = 90.0d0
    ENDIF
 
    ! CALL PHOTOLYSIS_KPP(zsd, I_ACF_in, I_ACF_out, TEMP_in, J_values)  ! Calculates J-values
    J_values = 1.0d-5
    J_values(1) = 6.0d-5  ! [s-1]

    IF (zsd_in >= 90) THEN
      J_values = 0.0d0
    ENDIF

    ! Some safety check
    DO i = 1,NPHOT
       IF (J_values(i) < -1d9 .OR. J_values(i) > 1d9 ) then
          WRITE(*,*) 'Note by Sampo: J_value(i) has bad value.'
          WRITE(*,*) 'i, J_value(i) = ', i, J_values(i)
          STOP
       ENDIF
    ENDDO

    CALL kpp_proceed(C_inout, TSTART_in, TEND_in, TEMP_in, RH, PRESS_in, O2, N2, M, H2O_in, RES1_in, RES2_in, J_values)

    RO2_out = RO2

  END SUBROUTINE chemistry

END MODULE chemistry_mod
