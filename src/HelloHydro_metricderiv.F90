#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating metric quantites 
subroutine HelloHydro_metricderiv(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer :: i,j,k
    cctk_real :: dtshift_term,dtbetacox,dtbetacoy,dtbetacoz
    
    call CCTK_INFO("Calculating derivatives of the 4-metric")
    
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                 
                ! derivs of tt component
                ! gttdt
                call calculate_dtshift_term(i,j,k,dtshift_term)
                gttdt(i,j,k) = -dtlapse(i,j,k) + dtshift_term
                ! No spatial variance so all other derivs are 0?
                gttdx(i,j,k) = 0.
                gttdy(i,j,k) = 0.
                gttdz(i,j,k) = 0. 
                ! derivs of tx component
                gtxdt(i,j,k) = 2.*dtbetacox  
                gtxdx(i,j,k) = 0.
                gtxdy(i,j,k) = 0.
                gtxdz(i,j,k) = 0.
                ! derivs of ty component 
                gtydt(i,j,k) = 2.*dtbetacoy  
                gtydx(i,j,k) = 0.
                gtydy(i,j,k) = 0.
                gtydz(i,j,k) = 0.
                ! derivs of tz component
                gtxdt(i,j,k) = 2.*dtbetacoz  
                gtxdx(i,j,k) = 0.
                gtxdy(i,j,k) = 0.
                gtxdz(i,j,k) = 0.

                ! derivs of xx component
                
            enddo 
        enddo 
    enddo
    
    print*, "Spatial derivs: ", gxxdx(1,1,1),gxxdy(1,1,1),gxxdz(1,1,1)
    
end subroutine HelloHydro_metricderiv

subroutine calculate_dtshift_term(i,j,k,dtshift_term)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in)    :: i,j,k
    cctk_real, intent(out) :: dtshift_term 
    cctk_real              :: betacox,betacoy,betacoz,dtbetacox,dtbetacoy,dtbetacoz

    ! Get the covariant component of the shfit vector
    call get_covshift(CCTK_ARGUMENTS,i,j,k,betacox,betacoy,betacoz)

    ! Get the covariant component of the derivatives of the shift vector
    call get_dtcovshift(CCTK_ARGUMENTS,i,j,k,dtbetacox,dtbetacoy,dtbetacoz)
    dtshift_term = betacox*dtbetax + betacoy*dtbetay + betacoz*dtbetaz &
     + betax*dtbetacox + betay*dtbetacoy + betaz*dtbetacoz
    
end subroutine calculate_dtshift_term