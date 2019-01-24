MODULE KPP_ROOT_Main
  ! A driver option for calling KPP as a module, from an external program
  ! by Sampo Smolander, 2008

  USE KPP_ROOT_Model
  USE KPP_ROOT_Initialize, ONLY: Initialize
  
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: KPP_SetUp, KPP_Proceed

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

  SUBROUTINE KPP_Proceed(CONS,TSTARTp,TENDp,TEMPp,RHp,PRESSp,O2p,N2p,Mp,H2Op,RES1p,RES2p,Jp,KVALUESp)

    IMPLICIT NONE
    REAL(kind=dp), INTENT(INOUT) :: CONS(NSPEC)
    REAL(kind=dp), INTENT(IN) :: TSTARTp,TENDp
    REAL(kind=dp), INTENT(IN) :: TEMPp,RHp,PRESSp,O2p,N2p,Mp,H2Op,RES1p,RES2p
    REAL(kind=dp), INTENT(IN) :: Jp(NPHOT),KVALUESp(42)

    REAL(kind=dp) :: T1,T2

!    T = TSTARTp
!    TIME = T
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

    KMT01   = KVALUESp(1)
    KMT02   = KVALUESp(2)
    KMT03   = KVALUESp(3)
    KMT04   = KVALUESp(4)
    KMT05   = KVALUESp(5)
    KMT06   = KVALUESp(6)
    KMT07   = KVALUESp(7)
    KMT08   = KVALUESp(8)
    KMT09   = KVALUESp(9)
    KMT10   = KVALUESp(10)
    KMT11   = KVALUESp(11)
    KMT12   = KVALUESp(12)
    KMT13   = KVALUESp(13)
    KMT14   = KVALUESp(14)
    KMT15   = KVALUESp(15)
    KMT16   = KVALUESp(16)
    KMT17   = KVALUESp(17)
    KMT18   = KVALUESp(18)
    KFPAN   = KVALUESp(19)
    KBPAN   = KVALUESp(20)
    KRO2NO  = KVALUESp(21)
    KRO2HO2 = KVALUESp(22)
    KAPHO2  = KVALUESp(23)
    KAPNO   = KVALUESp(24)
    KRO2NO3 = KVALUESp(25)
    KNO3AL  = KVALUESp(26)
    KDEC    = KVALUESp(27)
    KROPRIM = KVALUESp(28)
    KROSEC  = KVALUESp(29)
    KCH3O2  = KVALUESp(30)
    K298CH3O2 = KVALUESp(31)
    K_h2so4 = KVALUESp(32)
    K_CI = KVALUESp(33)
    KMT34 = KVALUESp(34)
    KMT35 = KVALUESp(35)
    KMT36 = KVALUESp(36)
    KMT37 = KVALUESp(37)
    KMT38 = KVALUESp(38)
    KMT39 = KVALUESp(39)
    KMT40 = KVALUESp(40)
    KMT41 = KVALUESp(41)
    KMT42 = KVALUESp(42)

    C = CONS
    ! Update_RCONST and INTEGRATE will operate on C
    CALL Update_RCONST()
    CALL INTEGRATE( TIN = T1, TOUT = T2, RSTATUS_U = RSTATE, ICNTRL_U = ICNTRL )
    CONS = C

!    T = RSTATE(1)
!    TIME = T
  END SUBROUTINE KPP_Proceed

END MODULE KPP_ROOT_Main
