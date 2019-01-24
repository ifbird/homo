PROGRAM main
  USE mcm_Parameters  ! NSPEC, ind_XXX, dp
  USE mcm_Monitor, ONLY: SPC_NAMES, EQN_NAMES
  USE chemistry_mod  ! chemistry, kpp_setup, NPHOT

  IMPLICIT NONE

  ! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
  REAL(dp), PARAMETER :: pi = 2*ASIN(1.)

  REAL(dp) :: c(NSPEC)
  REAL(dp) :: k
  REAL(dp) :: T, RH, PRES
  REAL(dp) :: O2, N2, M, H2O, H2
  REAL(dp) :: RO2
  REAL(dp) :: J_values(NPHOT)
  REAL(dp) :: I_ACF_in(75), I_ACF_out(75)

  REAL(dp) :: Ag, RGlob, zsd

  !***** time variables *****!
  REAL(dp) :: time, time_start, time_end
  REAL(dp) :: dt                          ! time step, in seconds
  REAL(dp) :: output_dt
  INTEGER  :: time_in_ms                  ! time in milliseconds
  INTEGER, PARAMETER :: one_hour = 60*60  ! in seconds

  ! Other parameters
  REAL(dp), PARAMETER :: ppm = 1.0d-6   ! [mol mol-1]
  REAL(dp), PARAMETER :: ppb = 1.0d-9   ! [mol mol-1]
  REAL(dp), PARAMETER :: ppt = 1.0d-12  ! [mol mol-1]
  REAL(dp), PARAMETER :: T00 = 273.15d0  ! [K]
  REAL(dp), PARAMETER :: P00 = 1.01325d5  ! [Pa]
  REAL(dp), PARAMETER :: Rgas = 8.314d0  ! [J mol-1 K-1]
  REAL(dp), PARAMETER :: NA = 6.022d23  ! [molec mol-1]

  ! Other variables
  REAL(dp) :: Nair = 2.25d19  ! [molec cm-3]


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Start program
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! write(*,*) NSPEC
  ! CALL InitSaveData()
  ! CALL Initialize()

  !------------------!
  ! Initial state
  !------------------!

  !!!!! Meteorology conditions
  T = -10.0d0 + T00  ! [K]
  RH = 0.70d0
  PRES = P00

  RGlob = 0.0d0
  zsd = 1.0d0
  Ag = 0.1d0
  I_ACF_in = 0.0d0

  !!!!! atmospheric O2, N2, H2O are kept constant
  O2 = 0.21d0*Nair  ! [molec cm-3]
  N2 = 0.78d0*Nair  ! [molec cm-3]
  M = N2 + O2  ! [molec cm-3]
  H2O = svp(T-T00)*RH*NA/(Rgas*T)*1.0d-6  ! [molec cm-3]

  !!!!! Time variables
  time_start = 0.0d0
  time_end = 960.0d0
  time = time_start
  dt = 0.5d0
  output_dt = 10.0d0

  ! Initial concentration
  c = 0.0d0
  ! c(ind_O1D) = 1d6  ! [molec cm-3], oxygen radical O concentration
  ! c(ind_O) = 1d6  ! [molec cm-3], oxygen radical O concentration
  c(ind_O3 ) = 122.0d0 *ppb*Nair  ! [molec cm-3]
  c(ind_DMS) = 3706.0d0*ppt*Nair  ! [molec cm-3]
  c(ind_NO2) = 10.0d0*ppt*Nair

  !!!!! Initialize KPP
  CALL kpp_setup()

  !!!!! Open a data file and write down the initial values
  WRITE(*,'(A8,F10.3,A6)'), 'time = ', time, '  seconds'
  OPEN(11,file="output/c.dat",status='replace',action='write')
  WRITE(11,*) c

  !------------------!
  ! Start time loop
  !------------------!
  DO WHILE (time < time_end)

    CALL chemistry(c, time, time+dt, T, RH, PRES, RGlob, zsd, Ag, &
      O2, N2, M, H2O, 1.0d-3, 1.0d-3, I_ACF_in, RO2, J_values, I_ACF_out)

    time = time + dt

    ! Write output every output_dt
    ! time_in_ms = FLOOR(1000*time)
    IF ( MODULO(FLOOR(time), FLOOR(output_dt)) == 0 ) THEN
      WRITE(*,'(A8,F10.3,A9)') 'time = ', time, '  seconds'
      WRITE(*,*) RO2
      WRITE(11,*) c
    END IF
    
  END DO

  CLOSE(11)


CONTAINS


  FUNCTION svp(TC)
    REAL(dp) :: TC  ! [degC]
    REAL(dp) :: svp  ! [Pa]

    REAL(dp), PARAMETER :: a0 = 6.107799961,     &
                           a1 = 4.436518524D-1,  &
                           a2 = 1.428945805D-2,  &
                           a3 = 2.650648471D-4,  &
                           a4 = 3.031240396D-6,  &
                           a5 = 2.034080948D-8,  &
                           a6 = 6.136820929D-11
    svp = (a0 + a1*TC + a2*TC**2 + a3*TC**3 + a4*TC**4 + a5*TC**5 + a6*TC**6) * 1.0d2
  END FUNCTION svp
END PROGRAM main
