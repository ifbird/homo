! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! The Stoichiometric Chemical Model File
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
! File                 : chapman_Stoichiom.f90
! Time                 : Sun Dec 11 22:03:00 2016
! Working directory    : /home/pzzhou/Documents/Research/Paper/Myself/2016/Iodine/kpp_iodine/Chemistry/test
! Equation file        : chapman.kpp
! Output root filename : chapman
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE chapman_Stoichiom

  USE chapman_Parameters
  USE chapman_StoichiomSP

  IMPLICIT NONE

CONTAINS


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! ReactantProd - Reactant Products in each equation
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      ARP       - Reactant product in each equation
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE ReactantProd ( V, F, ARP )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! ARP - Reactant product in each equation
  REAL(kind=dp) :: ARP(NREACT)


! Reactant Products in each equation are useful in the
!     stoichiometric formulation of mass action law
  ARP(1) = F(2)
  ARP(2) = V(2)*F(2)
  ARP(3) = V(3)
  ARP(4) = V(2)*V(3)
  ARP(5) = V(3)
  ARP(6) = V(1)*F(1)
  ARP(7) = V(1)*V(3)
  ARP(8) = V(3)*V(4)
  ARP(9) = V(2)*V(5)
  ARP(10) = V(5)
      
END SUBROUTINE ReactantProd

! End of ReactantProd function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! JacReactantProd - Jacobian of Reactant Products vector
!   Arguments :
!      V         - Concentrations of variable species (local)
!      F         - Concentrations of fixed species (local)
!      JVRP      - d ARP(1:NREACT)/d VAR (1:NVAR)
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

SUBROUTINE JacReactantProd ( V, F, JVRP )

! V - Concentrations of variable species (local)
  REAL(kind=dp) :: V(NVAR)
! F - Concentrations of fixed species (local)
  REAL(kind=dp) :: F(NFIX)
! JVRP - d ARP(1:NREACT)/d VAR (1:NVAR)
  REAL(kind=dp) :: JVRP(NJVRP)


! Reactant Products in each equation are useful in the
!    stoichiometric formulation of mass action law
! Below we compute the Jacobian of the Reactant Products vector
!    w.r.t. variable species: d ARP(1:NREACT) / d Var(1:NVAR)

! JVRP(1) = dARP(2)/dV(2)
  JVRP(1) = F(2)
! JVRP(2) = dARP(3)/dV(3)
  JVRP(2) = 1
! JVRP(3) = dARP(4)/dV(2)
  JVRP(3) = V(3)
! JVRP(4) = dARP(4)/dV(3)
  JVRP(4) = V(2)
! JVRP(5) = dARP(5)/dV(3)
  JVRP(5) = 1
! JVRP(6) = dARP(6)/dV(1)
  JVRP(6) = F(1)
! JVRP(7) = dARP(7)/dV(1)
  JVRP(7) = V(3)
! JVRP(8) = dARP(7)/dV(3)
  JVRP(8) = V(1)
! JVRP(9) = dARP(8)/dV(3)
  JVRP(9) = V(4)
! JVRP(10) = dARP(8)/dV(4)
  JVRP(10) = V(3)
! JVRP(11) = dARP(9)/dV(2)
  JVRP(11) = V(5)
! JVRP(12) = dARP(9)/dV(5)
  JVRP(12) = V(2)
! JVRP(13) = dARP(10)/dV(5)
  JVRP(13) = 1
      
END SUBROUTINE JacReactantProd

! End of JacReactantProd function
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



! Begin Derivative w.r.t. Rate Coefficients

! ------------------------------------------------------------------------------
! Subroutine for the derivative of Fun with respect to rate coefficients
! -----------------------------------------------------------------------------

      SUBROUTINE  dFun_dRcoeff( V, F, NCOEFF, JCOEFF, DFDR )
       
      USE chapman_Parameters
      USE chapman_StoichiomSP
      IMPLICIT NONE 

! V - Concentrations of variable/radical/fixed species            
      REAL(kind=dp) V(NVAR), F(NFIX)
! NCOEFF - the number of rate coefficients with respect to which we differentiate
      INTEGER NCOEFF       
! JCOEFF - a vector of integers containing the indices of reactions (rate
!          coefficients) with respect to which we differentiate
      INTEGER JCOEFF(NCOEFF)       
! DFDR  - a matrix containg derivative values; specifically, 
!         column j contains d Fun(1:NVAR) / d RCT( JCOEFF(j) )
!         for each 1 <= j <= NCOEFF
!         This matrix is stored in a column-wise linearized format
      REAL(kind=dp) DFDR(NVAR*NCOEFF)

! Local vector with reactant products
      REAL(kind=dp) A_RPROD(NREACT)
      REAL(kind=dp) aj
      INTEGER i,j,k
      
! Compute the reactant products of all reactions     
      CALL ReactantProd ( V, F, A_RPROD )

! Compute the derivatives by multiplying column JCOEFF(j) of the stoichiometric matrix with A_RPROD       
      DO j=1,NCOEFF
!                  Initialize the j-th column of derivative matrix to zero       
         DO i=1,NVAR
           DFDR(i+NVAR*(j-1)) = 0.0_dp 
         END DO
!                  Column JCOEFF(j) in the stoichiometric matrix times the
!                  reactant product  of the JCOEFF(j)-th reaction      
!                  give the j-th column of the derivative matrix   
         aj = A_RPROD(JCOEFF(j))
         DO k=CCOL_STOICM(JCOEFF(j)),CCOL_STOICM(JCOEFF(j)+1)-1
           DFDR(IROW_STOICM(k)+NVAR*(j-1)) = STOICM(k)*aj
         END DO
      END DO
      
      END SUBROUTINE  dFun_dRcoeff

! End Derivative w.r.t. Rate Coefficients


! Begin Jacobian Derivative w.r.t. Rate Coefficients

! ------------------------------------------------------------------------------
! Subroutine for the derivative of Jac with respect to rate coefficients
! Times a user vector
! -----------------------------------------------------------------------------

      SUBROUTINE  dJac_dRcoeff( V, F, U, NCOEFF, JCOEFF, DJDR )
       
      USE chapman_Parameters
      USE chapman_StoichiomSP
      IMPLICIT NONE 

! V - Concentrations of variable/fixed species            
      REAL(kind=dp) V(NVAR), F(NFIX)
! U - User-supplied Vector           
      REAL(kind=dp) U(NVAR)
! NCOEFF - the number of rate coefficients with respect to which we differentiate
      INTEGER NCOEFF       
! JCOEFF - a vector of integers containing the indices of reactions (rate
!          coefficients) with respect to which we differentiate
      INTEGER JCOEFF(NCOEFF)       
! DFDR  - a matrix containg derivative values; specifically, 
!         column j contains d Jac(1:NVAR) / d RCT( JCOEFF(j) ) * U
!                     for each 1 <= j <= NCOEFF
!         This matrix is stored in a column-wise linearized format
      REAL(kind=dp) DJDR(NVAR*NCOEFF)

! Local vector for Jacobian of reactant products
      REAL(kind=dp) JV_RPROD(NJVRP)
      REAL(kind=dp) aj
      INTEGER i,j,k
      
! Compute the Jacobian of all reactant products   
      CALL JacReactantProd( V, F, JV_RPROD )

! Compute the derivatives by multiplying column JCOEFF(j) of the stoichiometric matrix with A_PROD       
      DO j=1,NCOEFF
!                  Initialize the j-th column of derivative matrix to zero       
         DO i=1,NVAR
           DJDR(i+NVAR*(j-1)) = 0.0_dp
         END DO
!                  Column JCOEFF(j) in the stoichiometric matrix times the
!                  ( Gradient of reactant product of the JCOEFF(j)-th reaction X user vector )    
!                  give the j-th column of the derivative matrix   
!
!          Row JCOEFF(j) of JV_RPROD times the user vector
         aj = 0.0_dp
         DO k=CROW_JVRP(JCOEFF(j)),CROW_JVRP(JCOEFF(j)+1)-1
             aj = aj + JV_RPROD(k)*U(ICOL_JVRP(k))
         END DO
!          Column JCOEFF(j) of Stoichiom. matrix times aj         
         DO k=CCOL_STOICM(JCOEFF(j)),CCOL_STOICM(JCOEFF(j)+1)-1
           DJDR(IROW_STOICM(k)+NVAR*(j-1)) = STOICM(k)*aj
         END DO
      END DO
      
      END SUBROUTINE  dJac_dRcoeff

! End Jacobian Derivative w.r.t. Rate Coefficients


END MODULE chapman_Stoichiom

