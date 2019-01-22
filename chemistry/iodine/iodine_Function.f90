! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The ODE Function of Chemical Model File
! 
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
! 
! File                 : iodine_Function.f90
! Time                 : Tue Apr 24 18:06:11 2018
! Working directory    : /home/local/pzzhou/Documents/Research/Paper/Myself/2016/Iodine/kpp_iodine/CHEMISTRY/iodine
! Equation file        : iodine.kpp
! Output root filename : iodine
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE iodine_Function

  USE iodine_Parameters
  IMPLICIT NONE

! A - Rate for each equation
  REAL(kind=dp) :: A(NREACT)

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Fun - time derivatives of variables - Agregate form
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      RCT       - Rate constants (local)
!      Vdot      - Time derivative of variable species concentrations
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE Fun ( V, F, RCT, Vdot )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! RCT - Rate constants (local)
  REAL(kind=dp) :: RCT(NREACT)
! Vdot - Time derivative of variable species concentrations
  REAL(kind=dp) :: Vdot(NVAR)


! Computation of equation rates
  A(1) = RCT(1)*V(23)
  A(2) = RCT(2)*V(14)
  A(3) = RCT(3)*V(21)
  A(4) = RCT(4)*V(21)*V(23)
  A(5) = RCT(5)*V(14)*V(21)
  A(6) = RCT(6)*V(14)*V(21)
  A(7) = RCT(7)*V(14)*V(19)
  A(8) = RCT(8)*V(26)*V(27)
  A(9) = RCT(9)*V(15)*V(26)
  A(10) = RCT(10)*V(21)*V(26)
  A(11) = RCT(11)*V(26)*V(26)
  A(12) = RCT(12)*V(10)
  A(13) = RCT(13)*V(23)*V(26)
  A(14) = RCT(14)*V(23)*V(27)
  A(15) = RCT(15)*V(26)*V(26)
  A(16) = RCT(16)*V(15)*V(23)
  A(17) = RCT(17)*V(15)*V(26)
  A(18) = RCT(18)*V(21)*V(26)
  A(19) = RCT(19)*V(27)*V(27)
  A(20) = RCT(20)*V(21)*V(27)
  A(21) = RCT(21)*V(23)*V(24)
  A(22) = RCT(22)*V(23)*V(29)
  A(23) = RCT(23)*V(23)*V(29)
  A(24) = RCT(24)*V(23)*V(25)
  A(25) = RCT(25)*V(11)*V(26)
  A(26) = RCT(26)*V(26)*V(29)
  A(27) = RCT(27)*V(26)*V(29)
  A(28) = RCT(28)*V(25)*V(26)
  A(29) = RCT(29)*V(24)*V(27)
  A(30) = RCT(30)*V(27)*V(29)
  A(31) = RCT(31)*V(25)*V(27)
  A(32) = RCT(32)*V(25)*V(27)
  A(33) = RCT(33)*V(21)*V(24)
  A(34) = RCT(34)*V(24)*V(29)
  A(35) = RCT(35)*V(24)*V(25)
  A(36) = RCT(36)*V(29)*V(29)
  A(37) = RCT(37)*V(25)*V(29)
  A(38) = RCT(38)*V(20)*V(23)
  A(39) = RCT(39)*V(23)*V(28)
  A(40) = RCT(40)*V(22)*V(27)
  A(41) = RCT(41)*V(21)*V(22)
  A(42) = RCT(42)*V(22)*V(24)
  A(43) = RCT(43)*V(22)*V(29)
  A(44) = RCT(44)*V(22)*V(25)
  A(45) = RCT(45)*V(20)*V(25)
  A(46) = RCT(46)*V(17)*V(26)
  A(47) = RCT(47)*V(20)*V(26)
  A(48) = RCT(48)*V(25)*V(28)
  A(49) = RCT(49)*V(27)*V(28)
  A(50) = RCT(50)*V(21)*V(28)
  A(51) = RCT(51)*V(21)*V(28)
  A(52) = RCT(52)*V(28)*V(28)
  A(53) = RCT(53)*V(28)*V(28)
  A(54) = RCT(54)*V(28)*V(28)
  A(55) = RCT(55)*V(28)*V(28)
  A(56) = RCT(56)*V(24)*V(28)
  A(57) = RCT(57)*V(28)*V(29)
  A(58) = RCT(58)*V(12)*V(12)
  A(59) = RCT(59)*V(13)*V(13)
  A(60) = RCT(60)*V(18)*V(24)
  A(61) = RCT(61)*V(17)*V(28)
  A(62) = RCT(62)*V(16)*V(26)
  A(63) = RCT(63)*V(18)*V(26)
  A(64) = RCT(64)*V(27)*V(28)

! Aggregate function
  Vdot(1) = A(4)+A(6)+A(13)+A(18)+A(27)+A(31)+A(39)+A(40)+A(41)+A(49)+A(50)
  Vdot(2) = A(26)
  Vdot(3) = A(30)
  Vdot(4) = A(34)
  Vdot(5) = A(36)
  Vdot(6) = A(37)
  Vdot(7) = A(45)+A(57)
  Vdot(8) = A(55)
  Vdot(9) = A(63)+A(64)
  Vdot(10) = -A(12)+A(13)
  Vdot(11) = -A(25)
  Vdot(12) = A(42)-2*A(58)
  Vdot(13) = A(43)-2*A(59)
  Vdot(14) = -A(2)+A(3)-A(5)-A(6)-A(7)
  Vdot(15) = -A(9)+A(15)-A(16)-A(17)+A(19)
  Vdot(16) = A(47)+A(49)+A(61)-A(62)
  Vdot(17) = A(40)-A(46)-A(61)
  Vdot(18) = A(48)+A(51)+A(54)-A(60)-A(63)
  Vdot(19) = -A(7)+A(8)+A(11)+A(17)+A(25)+A(46)+A(62)
  Vdot(20) = -A(38)-A(45)-A(47)+A(52)+A(58)+A(59)
  Vdot(21) = A(1)-A(3)-A(4)-A(5)-A(6)-A(10)-A(18)-A(20)-A(33)-A(41)-A(50)-A(51)
  Vdot(22) = A(38)+A(39)-A(40)-A(41)-A(42)-A(43)-A(44)+A(45)+A(46)+A(47)+A(50)+2*A(53)+A(54)+A(56)+A(61)
  Vdot(23) = -A(1)+A(2)-A(4)+2*A(5)+A(11)-A(13)-A(14)-A(16)-A(21)-A(22)-A(23)-A(24)-A(38)-A(39)
  Vdot(24) = -A(21)+A(22)-A(29)-A(33)-A(34)-A(35)-A(42)-A(56)+2*A(58)-A(60)
  Vdot(25) = A(23)-A(24)-A(28)-A(31)-A(32)-A(35)-A(37)-A(44)-A(45)-A(48)
  Vdot(26) = 2*A(7)-A(8)-A(9)-A(10)-2*A(11)-A(13)+A(14)-2*A(15)+A(16)-A(17)-A(18)+A(20)-A(25)-A(26)-A(27)-A(28)+A(29)&
               &+A(32)-A(46)-A(47)-A(62)-A(63)
  Vdot(27) = -A(8)+A(9)+A(10)+A(12)-A(14)+A(16)+A(17)+A(18)-2*A(19)-A(20)+A(28)-A(29)-A(30)-A(31)-A(32)-A(40)-A(49)&
               &-A(64)
  Vdot(28) = A(38)-A(39)+A(41)+A(44)-A(48)-A(49)-A(50)-A(51)-2*A(52)-2*A(53)-2*A(54)-2*A(55)-A(56)-A(57)+A(60)-A(61)&
               &+A(62)-A(64)
  Vdot(29) = A(21)-A(22)-A(23)+A(24)+A(25)-A(26)-A(27)+A(28)+A(29)-A(30)+A(32)+A(33)-A(34)+2*A(35)-2*A(36)-A(37)-A(43)&
               &+A(44)+A(48)+A(56)-A(57)+2*A(59)+A(60)
      
END SUBROUTINE Fun

! End of Fun function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



END MODULE iodine_Function
