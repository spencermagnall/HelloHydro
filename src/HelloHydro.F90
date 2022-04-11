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
    CCTK_REAL :: PI, rho1, rho2, rho3
    integer   :: i,j,k
    logical :: INIT = .TRUE. 

    !if (INIT .eqv. .TRUE.) then 
    !    rho = 13.29
    !    INIT = .FALSE.
    !endif
    PI = 4.D0*DATAN(1.D0)
    print*, rho(1,1,1)
    rho1 = rho(1,1,1)
    rho2 = rho_p(1,1,1)
    rho3 = rho_p_p(1,1,1)
    !print*, "eTtt = ", eTtt
    !allocate(hub(size(rho),size(rho),size(rho)))
    ! Hubble parameter, comes from friedman     equation for dust and k = 0  
    print*, "rho values are: ", rho1,rho2,rho3
    ! Get the hubble parameter
    call get_hub(CCTK_ARGUMENTS,hub) 
    
    print *, "Hub valvue is: ", hub
    !call CCTK_OutputVar(rho1)
    call CCTK_INFO("Running Hello Hydro")
    !CCTK_REAL :: rho,press,eps,vel[3]
    ! One liner version of Phantom
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)

                    ! rho_p is rho at t rho is rho at t+dt 
                    rho(i,j,k) = rho_p(i,j,k) - 3*hub*cctk_delta_time*rho_p(i,j,k)


            enddo 
        enddo 
    enddo 
    !rho = 13.29
    !call CCTK_INFO(rho)
    !Phantom stuff here 
    ! Call evol_step
    ! What needs to go into evol_step 
    ! time so that phantom knows when to dump
    ! dt so that particles can be stepped properly
    ! Which dt do we use? i.e does it need to be calculated in a seperate    phantom subroutine and passed to ET
    ! Do we need dtext?, dtnew? What do this variables even do?
    ! dtnew is the new timestep for the next itteration
end subroutine HelloHydro

subroutine get_hub(CCTK_ARGUMENTS,hub)
    implicit none 
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    CCTK_REAL, intent(out) :: hub
    CCTK_REAL :: PI
    CCTK_REAL :: rho1

    PI = 4.D0*DATAN(1.D0)
    ! FIX THIS UP 
    if (timelevels .eq. 1) then 
        rho1 = rho(1,1,1)
    elseif(timelevels .lt. 2) then 
        rho1 = rho_p(1,1,1)
    else
        rho1 = rho_p_p(1,1,1)
    endif 
    hub =  sqrt((8.*PI*rho1/3.))
end subroutine 
