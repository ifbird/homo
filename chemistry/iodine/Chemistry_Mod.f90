! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Chemistry module, created from ozone_Main.f90.
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

MODULE Chemistry_Mod

USE iodine_Model
USE iodine_Initialize, ONLY: Initialize


IMPLICIT NONE

PRIVATE
PUBLIC :: kpp_initialize, kpp_chemistry

REAL(kind=dp) :: RSTATE(20)
INTEGER :: I


CONTAINS


SUBROUTINE kpp_initialize()

  STEPMIN = 0.0d0
  STEPMAX = 0.0d0

  DO I=1,NVAR
    RTOL(I) = 1.0d-4
    ATOL(I) = 1.0d-3
  END DO

  CALL Initialize()
END SUBROUTINE kpp_initialize


SUBROUTINE kpp_chemistry(C_in, TSTART_in, TEND_in, TEMP_in, RH_in, PRESS_in, O2_in, N2_in, k_out)
  REAL(dp), INTENT(INOUT) :: C_in(NSPEC)  ! NSPEC: total number of species, NVAR: number of variable species
  REAL(dp), INTENT(IN) :: TSTART_in, TEND_in
  REAL(dp), INTENT(IN) :: TEMP_in, RH_in, PRESS_in, O2_in, N2_in

  REAL(dp), INTENT(OUT) :: k_out

  REAL(dp) :: T1, T2

  T1 = TSTART_in
  T2 = TEND_in
  TIME = T1  ! current time, used globally for KPP modules

  TEMP = TEMP_in
  O2 = O2_in
  N2 = N2_in

  C(1:NVAR) = C_in(1:NVAR)  ! Do not set fixed variables from outside

  CALL calculate_k()
  CALL Update_SUN()
  CALL Update_RCONST()

  CALL INTEGRATE( TIN = T1, TOUT = T2, RSTATUS_U = RSTATE, &
    ICNTRL_U = (/ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 /) )

  C_in = C

END SUBROUTINE kpp_chemistry

SUBROUTINE calculate_k()
  R01 = 0.00831d0  ! kJ K-1 mol-1
  K01 = 0.0d0
END SUBROUTINE calculate_k

! SUBROUTINE calculate_k(time)
!   ! This subroutine will write into the global variable 'kval'
!   REAL(dp), INTENT(in) :: time
!   REAL(dp), PARAMETER  :: c3 = 22.62d0
!   REAL(dp), PARAMETER  :: c4 = 7.601d0
!   REAL(dp), PARAMETER  :: pi = 2*ASIN(1.0d0)
!   REAL(dp), PARAMETER  :: w = pi/(12*3600)
!   REAL(dp)             :: coswt
! 
!   coswt = COS(w*time)
!   IF (coswt < 0) THEN
!     k(1) = EXP(c3/coswt)
!     k(2) = EXP(c4/coswt)
!   ELSE
!     k(1) = 0.0d0
!     k(2) = 0.0d0
!   END IF
! END SUBROUTINE calculate_k


END MODULE Chemistry_Mod

! End of MAIN function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
