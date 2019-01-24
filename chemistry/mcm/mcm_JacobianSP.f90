! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! 
! Sparse Jacobian Data Structures File
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
! File                 : mcm_JacobianSP.f90
! Time                 : Thu Jan 24 12:11:47 2019
! Working directory    : /home/local/pzzhou/Scripts/repos/homo/chemistry/mcm
! Equation file        : mcm.kpp
! Output root filename : mcm
! 
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



MODULE mcm_JacobianSP

  PUBLIC
  SAVE


! Sparse Jacobian Data


  INTEGER, PARAMETER, DIMENSION(360) :: LU_IROW_0 = (/ &
       1,  1,  1,  1,  1,  2,  2,  2,  3,  3,  4,  4, &
       4,  4,  4,  5,  5,  6,  6,  6,  7,  7,  7,  8, &
       8,  8,  9,  9,  9, 10, 10, 10, 11, 11, 11, 11, &
      12, 12, 12, 13, 13, 13, 14, 14, 14, 15, 15, 15, &
      16, 16, 16, 17, 17, 17, 17, 17, 18, 18, 18, 18, &
      18, 19, 19, 19, 19, 20, 20, 20, 20, 21, 21, 21, &
      21, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, &
      24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, &
      26, 26, 26, 26, 26, 26, 27, 27, 27, 28, 28, 28, &
      28, 29, 29, 29, 29, 30, 30, 30, 30, 31, 31, 31, &
      31, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33, &
      33, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, &
      35, 36, 36, 36, 36, 36, 37, 37, 37, 37, 38, 38, &
      38, 38, 38, 38, 38, 38, 38, 38, 39, 39, 39, 39, &
      39, 39, 39, 40, 40, 40, 40, 40, 40, 40, 40, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, &
      41, 41, 41, 42, 42, 42, 42, 42, 42, 43, 43, 43, &
      43, 43, 43, 43, 44, 44, 44, 44, 44, 44, 44, 44, &
      44, 44, 44, 44, 45, 45, 45, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 47, 47, 47, 47, &
      47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 47, 48, &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, 48, &
      48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
      50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52 /)
  INTEGER, PARAMETER, DIMENSION(81) :: LU_IROW_1 = (/ &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 53, 53, 53, 53, 53, 53, &
      53, 53, 53, 53, 53, 53, 53, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55 /)
  INTEGER, PARAMETER, DIMENSION(441) :: LU_IROW = (/&
    LU_IROW_0, LU_IROW_1 /)

  INTEGER, PARAMETER, DIMENSION(360) :: LU_ICOL_0 = (/ &
       1, 39, 49, 52, 53,  2,  7, 33,  3,  4,  4,  6, &
      38, 39, 47,  5, 53,  6, 47, 49,  7, 50, 54,  8, &
      49, 52,  9, 41, 49, 10, 45, 54, 11, 36, 49, 51, &
      12, 45, 49, 13, 42, 49, 14, 49, 51, 15, 36, 49, &
      16, 40, 49, 17, 31, 42, 50, 51, 18, 38, 46, 49, &
      52, 19, 49, 52, 55, 20, 49, 52, 54, 21, 46, 49, &
      52, 22, 23, 24, 41, 49, 50, 16, 23, 34, 40, 49, &
      13, 24, 31, 42, 49, 16, 25, 34, 40, 49, 50, 51, &
      26, 30, 32, 45, 50, 51, 27, 49, 50, 28, 46, 49, &
      54, 29, 49, 54, 55, 30, 45, 49, 51, 31, 42, 49, &
      52, 32, 45, 49, 52, 27, 28, 33, 41, 46, 49, 50, &
      54, 34, 40, 49, 52, 15, 29, 35, 36, 49, 52, 54, &
      55, 27, 36, 49, 50, 51, 37, 44, 51, 54, 18, 21, &
      38, 46, 48, 49, 50, 51, 52, 53,  5, 39, 47, 50, &
      51, 53, 54, 11, 34, 36, 40, 49, 50, 51, 52, 12, &
      17, 25, 26, 30, 31, 32, 34, 40, 41, 42, 45, 49, &
      50, 51, 52, 27, 42, 49, 50, 51, 52, 37, 43, 44, &
      51, 53, 54, 55, 17, 24, 31, 37, 42, 44, 49, 50, &
      51, 52, 53, 54, 10, 15, 32, 35, 36, 37, 38, 43, &
      44, 45, 46, 48, 49, 50, 51, 52, 53, 54, 55, 21, &
      28, 46, 48, 49, 50, 51, 52, 54, 35, 36, 37, 39, &
      43, 44, 47, 48, 49, 50, 51, 52, 53, 54, 55, 19, &
      23, 25, 34, 40, 43, 44, 46, 48, 49, 50, 51, 52, &
      53, 54, 55,  5,  8,  9, 11, 12, 13, 14, 15, 16, &
      18, 19, 20, 21, 22, 23, 24, 27, 28, 29, 30, 31, &
      32, 33, 34, 35, 36, 38, 40, 41, 42, 45, 46, 47, &
      48, 49, 50, 51, 52, 53, 54, 55,  7, 27, 33, 37, &
      39, 40, 41, 42, 44, 45, 46, 47, 48, 49, 50, 51, &
      52, 53, 54, 55, 14, 36, 37, 39, 40, 42, 43, 44, &
      45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55,  6 /)
  INTEGER, PARAMETER, DIMENSION(81) :: LU_ICOL_1 = (/ &
       8,  9, 12, 13, 16, 20, 22, 23, 24, 26, 30, 31, &
      32, 34, 36, 38, 40, 41, 42, 45, 46, 47, 48, 49, &
      50, 51, 52, 53, 54, 55, 39, 43, 44, 46, 47, 48, &
      49, 50, 51, 52, 53, 54, 55,  7, 10, 14, 20, 28, &
      29, 30, 33, 36, 37, 39, 40, 41, 42, 43, 44, 45, &
      46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 19, 29, &
      43, 44, 49, 50, 51, 52, 53, 54, 55 /)
  INTEGER, PARAMETER, DIMENSION(441) :: LU_ICOL = (/&
    LU_ICOL_0, LU_ICOL_1 /)

  INTEGER, PARAMETER, DIMENSION(56) :: LU_CROW = (/ &
       1,  6,  9, 11, 16, 18, 21, 24, 27, 30, 33, 37, &
      40, 43, 46, 49, 52, 57, 62, 66, 70, 74, 80, 85, &
      90, 97,103,106,110,114,118,122,126,134,138,146, &
     151,155,165,172,180,196,202,209,221,240,249,264, &
     280,321,341,360,391,404,431,442 /)

  INTEGER, PARAMETER, DIMENSION(56) :: LU_DIAG = (/ &
       1,  6,  9, 11, 16, 18, 21, 24, 27, 30, 33, 37, &
      40, 43, 46, 49, 52, 57, 62, 66, 70, 74, 81, 86, &
      91, 97,103,106,110,114,118,122,128,134,140,147, &
     151,157,166,175,189,197,203,214,230,242,255,272, &
     314,335,355,387,401,429,441,442 /)


END MODULE mcm_JacobianSP

