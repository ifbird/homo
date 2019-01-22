PROGRAM main
  USE iodine_Parameters  ! NSPEC, ind_XXX, iodine_Precision
  USE iodine_Monitor  ! SPC_NAMES, EQN_NAMES
  USE Chemistry_Mod


  IMPLICIT NONE

  ! INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15,307)
  REAL(dp), PARAMETER :: pi = 2*ASIN(1.)

  REAL(dp) :: c(NSPEC)
  REAL(dp) :: k
  REAL(dp) :: T, RH, PRES
  REAL(dp) :: O2, N2, M

  !***** time variables *****!
  REAL(dp) :: time, time_start, time_end
  REAL(dp) :: dt                          ! time step, in seconds
  REAL(dp) :: output_dt
  INTEGER  :: time_in_ms                  ! time in milliseconds
  INTEGER, PARAMETER :: one_hour = 60*60  ! in seconds


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! Start program
  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
  ! write(*,*) NSPEC
  ! CALL InitSaveData()
  ! CALL Initialize()

  !------------------!
  ! Initial state
  !------------------!
  c = 0.0d0
  c(ind_O) = 1d6  ! [molec cm-3], oxygen radical O concentration
  c(ind_O3) = 1d12  ! [molec cm-3]
  c(ind_H2O) = 1.0d16  ! [molec cm-3]

  !!!!! atmospheric O2, N2, H2O are kept constant
  O2 = 0.21*2.4D19  ! [molec cm-3]
  N2 = 0.78*2.4D19  ! [molec cm-3]
  M = N2 + O2  ! [molec cm-3]

  !!!!! Meteorology conditions
  T = 290.0d0
  RH = 0.80d0
  PRES = 1.0d5

  !!!!! Time variables
  time_start = 12*3600.0d0
  time_end = 13*3600.0d0
  time = time_start
  dt = 1.0d0
  output_dt = 5.0d0

  !!!!! Initialize KPP
  CALL kpp_initialize()

  !!!!! Open a data file and write down the initial values
  WRITE(*,'(A8,F10.3,A6)'), 'time = ', time, '  seconds'
  OPEN(11,file="OUTPUT/c.dat",status='replace',action='write')
  WRITE(11,*) c

  !------------------!
  ! Start time loop
  !------------------!
  DO WHILE (time < time_end)

    CALL kpp_chemistry(c, time, time+dt, T, RH, PRES, O2, N2, k)

    time = time + dt

    ! Write output every output_dt
    ! time_in_ms = FLOOR(1000*time)
    IF ( MODULO(FLOOR(time), FLOOR(output_dt)) == 0 ) THEN  ! what a hack
      WRITE(*,'(A8,F10.3,A9)'), 'time = ', time, '  seconds'
      WRITE(11,*) c
    END IF
    
  END DO

  CLOSE(11)
END PROGRAM main
