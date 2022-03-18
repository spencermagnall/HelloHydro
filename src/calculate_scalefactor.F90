#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine calculate_scalefactor(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    REAL :: a,a3,a6

    call CCTK_INFO("Calculating Scale factor from metric determinant")
    a6 = -(gxx*gyy*gzz)
    a3 = a6**(1./2.)
    a = a3**(1./3.)
end subroutine calculate_scalefactor 
