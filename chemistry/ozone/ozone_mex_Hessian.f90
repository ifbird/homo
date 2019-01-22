
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 SUBROUTINE mexFunction( nlhs, plhs, nrhs, prhs )
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!                  Matlab Gateway for the Function Hessian
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 USE ozone_Model

      INTEGER nlhs, nrhs
      INTEGER plhs(*), prhs(*)
      INTEGER mxGetPr, mxCreateFull, mxGetM, mxgetN
      INTEGER VPtr, FPtr, RPtr, HESSPtr
      REAL(kind=dp) V(2), F(1), RCT(4)
      REAL(kind=dp) HESS(2)

! Check for the right number of input arguments
      IF ( nrhs .ne. 3 ) THEN
         CALL mexErrMsgTxt('Hessian requires 3 input vectors: &
     &V(2), F(1), RCT(4)')
      END IF 
! Check for the right number of output arguments
      IF ( nlhs .ne. 1 ) THEN
         CALL mexErrMsgTxt('Hessian requires 1 output vector: &
     &HESS(2)')
      END IF 

      plhs(1) = mxCreateDoubleMatrix(2,1,0)

      VPtr = mxGetPr(prhs(1));
      CALL mxCopyPtrToReal8(VPtr,V,2)
      
      FPtr = mxGetPr(prhs(2));
      CALL mxCopyPtrToReal8(FPtr,F,1)
      
      RPtr = mxGetPr(prhs(3));
      CALL mxCopyPtrToReal8(RPtr,RCT,4)

      HESSPtr = mxGetPr(plhs(1))

      CALL Hessian( V, F, RCT, HESS )

      CALL mxCopyReal8ToPtr(HESS, HESSPtr, 2)

 END SUBROUTINE mexFunction
