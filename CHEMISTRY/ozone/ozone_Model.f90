MODULE ozone_Model

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Completely defines the model ozone
!    by using all the associated modules
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  USE ozone_Precision
  USE ozone_Parameters
  USE ozone_Global
  USE ozone_Function
  USE ozone_Integrator
  USE ozone_Rates
  USE ozone_Jacobian
  USE ozone_Hessian
  USE ozone_Stoichiom
  USE ozone_LinearAlgebra
  USE ozone_Monitor
  USE ozone_Util

END MODULE ozone_Model

