#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine HelloHydro_init(CCTK_ARGUMENTS)
    implicit none
    integer :: i,j,k
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    call CCTK_INFO("Setting rho to be 0.")
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)

                rho(i,j,k) = 0.

            enddo 
        enddo 
    enddo 

end subroutine HelloHydro_init
