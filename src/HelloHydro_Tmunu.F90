#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating the stress energy tensor and adding its components

! Eventually this should be similar to GRHydro_Tmunu and calculate Tmunu properly but for now we just use dust solution

! I assume that everything is MCRF and that the 4 velocity is just u0 = 1, but this could be incorrect?? 
subroutine HelloHydro_Tmunu(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS

    call CCTK_INFO("Calculating components of stress energy tensor from primative variables")
    print*, "Value of Stress energy tensor before is: ", eTtt(1,1,1)
    eTtt = eTtt + rho
    eTxx = eTxx + press
    eTyy = eTyy + press
    eTzz = eTzz + press
    print*, "Value of Stress energy tensor is: ", eTtt(1,1,1)
    print*, "values of metric are: ", gxx(1,1,1), gyy(1,1,1), gzz(1,1,1)
    !stop 
end subroutine HelloHydro_Tmunu 
