#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating metric quantites 
subroutine HelloHydro_metricderiv(CCTK_ARGUMENTS)
    use metric_utils
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer :: i,j,k
    CCTK_REAL :: dtshift_term,dtbetacox,dtbetacoy,dtbetacoz
    CCTK_REAL :: liederiv(3) ! Not sure what the correct rank for this is 


    ! Instead of one big derivative array store time and spatial derivs
    ! as these are what is required for the source terms in phantom 

    
    call CCTK_INFO("Calculating derivatives of the 4-metric")

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Spatial derivatives 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Calculate the spatial derivs of the 3 metric with finite diff
    ! Should already by scheduled with central finite dont need to calculate here
    !call central_finite(CCTK_ARGUMENTS)
    !print*, "Finished finite differencing!"
    ! Set the tt and ti component spatial derivs to 0.
    ! tt comp
    gttdx = 0.
    gttdy = 0.
    gttdz = 0.
    ! tx comp
    gtxdx = 0.
    gtxdy = 0.
    gtxdz = 0.
    ! ty comp
    gtydx = 0.
    gtydy = 0.
    gtydz = 0.
    ! tz comp
    gtzdx = 0.
    gtzdy = 0.
    gtzdz = 0.
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Time derivatives 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    liederiv = 0.
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                 
                ! derivs of tt component
                ! gttdt
                call calculate_dtshift_term(CCTK_ARGUMENTS,i,j,k,dtshift_term)
                gttdt(i,j,k) = -dtalp(i,j,k) + dtshift_term
                
                ! Get the covariant component of the derivatives of the shift vector
                call get_dtcovshift(CCTK_ARGUMENTS,i,j,k,dtbetacox,dtbetacoy,dtbetacoz)

                ! derivs of tx component
                gtxdt(i,j,k) = 2.*dtbetacox  
               
                ! derivs of ty component 
                gtydt(i,j,k) = 2.*dtbetacoy
                
                ! derivs of tz component
                gtzdt(i,j,k) = 2.*dtbetacoz

                ! derivs of xx component
                gxxdt(i,j,k) = -2*alp(i,j,k)*kxx(i,j,k)

                ! derivs of xy component
                gxydt(i,j,k) = -2*alp(i,j,k)*kxy(i,j,k)

                ! derivs of xz component
                gxzdt(i,j,k) = -2*alp(i,j,k)*kxz(i,j,k)

                ! derivs of yy component
                gyydt(i,j,k) = -2*alp(i,j,k)*kyy(i,j,k)

                ! derivs of yz component
                gyzdt(i,j,k) = -2*alp(i,j,k)*kyz(i,j,k)

                 ! derivs of zz component
                gzzdt(i,j,k) = -2*alp(i,j,k)*kzz(i,j,k)

            enddo 
        enddo 
    enddo
    print*, "Metric time derivatives are : ", gttdt(2,2,2), gtxdt(2,2,2), gtydt(2,2,2), gtzdt(2,2,2), gxxdt(2,2,2),gxydt(2,2,2),gxzdt(2,2,2)
    print*, gyydt(2,2,2), gyzdt(2,2,2), gzzdt(2,2,2)
    stop
end subroutine HelloHydro_metricderiv

subroutine calculate_dtshift_term(CCTK_ARGUMENTS,i,j,k,dtshift_term)
    use metric_utils
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    
    integer, intent(in)    :: i,j,k
    CCTK_REAL, intent(out) :: dtshift_term 
    CCTK_REAL             :: betacox,betacoy,betacoz,dtbetacox,dtbetacoy,dtbetacoz

    ! Get the covariant component of the shfit vector
    call get_covshift(CCTK_ARGUMENTS,i,j,k,betacox,betacoy,betacoz)

    ! Get the covariant component of the derivatives of the shift vector
    call get_dtcovshift(CCTK_ARGUMENTS,i,j,k,dtbetacox,dtbetacoy,dtbetacoz)
    dtshift_term = betacox*dtbetax(i,j,k) + betacoy*dtbetay(i,j,k) + betacoz*dtbetaz(i,j,k) &
     + betax(i,j,k)*dtbetacox + betay(i,j,k)*dtbetacoy + betaz(i,j,k)*dtbetacoz
    
end subroutine calculate_dtshift_term