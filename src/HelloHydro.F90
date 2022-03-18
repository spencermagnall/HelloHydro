#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine HelloHydro(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    !CCTK_REAL :: INTENT(INOUT) :: rho
    CCTK_REAL :: hub  
    CCTK_REAL :: PI, rho1
    integer   :: i,j,k
    PI = 4.D0*DATAN(1.D0)
    print*, rho(1,1,1)
    rho1 = rho(1,1,1)
    !print*, "eTtt = ", eTtt
    !allocate(hub(size(rho),size(rho),size(rho)))
    ! Hubble parameter, comes from friedman     equation for dust and k = 0  
    print*, "rho1 value is: ", rho1
    hub = sqrt((8.*PI*rho1/3.))
    
    print *, "Hub valvue is: ", hub
    !call CCTK_OutputVar(rho1)
    call CCTK_INFO("Running Hello Hydro")
    !CCTK_REAL :: rho,press,eps,vel[3]
    ! One liner version of Phantom
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)

                rho(i,j,k) = rho(i,j,k) - 3*hub*cctk_delta_time*rho(i,j,k)
            enddo 
        enddo 
    enddo 
    !rho = 13.29
    !call CCTK_INFO(rho)
end subroutine HelloHydro 
