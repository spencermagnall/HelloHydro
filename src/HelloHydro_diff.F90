#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating metric quantites 
subroutine HelloHydro_diff(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer :: i,j,k
    
    call CCTK_INFO("Calculating derivatives of the spatial metric using finite difference")
    
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                 
                ! Calculate spatial derivatives using finite difference 
                call central_finite(CCTK_ARGUMENTS,i,j,k)
                
            enddo 
        enddo 
    enddo
    
    print*, "Spatial derivs: ", gxxdx(1,1,1),gxxdy(1,1,1),gxxdz(1,1,1)
    
end subroutine HelloHydro_diff

subroutine central_finite(CCTK_ARGUMENTS,i,j,k)
    implicit none 
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: i,j,k

    
    ! Get first derivatives of the spatial metric using a centered finite difference 
    ! Refactor later into a new function as there is a lot of copy and paste here 
    ! Might also be a better idea to use a general difff class to alow the code to be more modular 
    
    ! Periodic boundary conditions are not handled, so the derivatives on the boundaries are
    ! Incorrect

    ! gxx 
    gxxdx(i,j,k) = (gxx(i+1,j,k) - gxx(i-1,j,k))/(2.*cctk_delta_space(1))
    print*, "gxx values: ", gxx(i+1,j,k), gxx(i-1,j,k), gxxdx(i,j,k)
    stop
    gxxdy(i,j,k) = (gxx(i,j+1,k) - gxx(i,j-1,k))/(2.*cctk_delta_space(2))
    gxxdz(i,j,k) = (gxx(i,j,k+1) - gxx(i,j,k-1))/(2.*cctk_delta_space(3))
    ! gxy 
    gxydx(i,j,k) = (gxy(i+1,j,k) - gxy(i-1,j,k))/(2.*cctk_delta_space(1))
    gxydy(i,j,k) = (gxy(i,j+1,k) - gxy(i,j-1,k))/(2.*cctk_delta_space(2))
    gxydz(i,j,k) = (gxy(i,j,k+1) - gxy(i,j,k-1))/(2.*cctk_delta_space(3))
    ! gxz
    gxzdx(i,j,k) = (gxz(i+1,j,k) - gxz(i-1,j,k))/(2.*cctk_delta_space(1))
    gxzdy(i,j,k) = (gxz(i,j+1,k) - gxz(i,j-1,k))/(2.*cctk_delta_space(2))
    gxzdz(i,j,k) = (gxz(i,j,k+1) - gxz(i,j,k-1))/(2.*cctk_delta_space(3))
    ! gyy 
    gyydx(i,j,k) = (gyy(i+1,j,k) - gyy(i-1,j,k))/(2.*cctk_delta_space(1))
    gyydy(i,j,k) = (gyy(i,j+1,k) - gyy(i,j-1,k))/(2.*cctk_delta_space(2))
    gyydz(i,j,k) = (gyy(i,j,k+1) - gyy(i,j,k-1))/(2.*cctk_delta_space(3))

    ! gyz 
    gyzdx(i,j,k) = (gyz(i+1,j,k) - gyz(i-1,j,k))/(2.*cctk_delta_space(1))
    gyzdy(i,j,k) = (gyz(i,j+1,k) - gyz(i,j-1,k))/(2.*cctk_delta_space(2))
    gyzdz(i,j,k) = (gyz(i,j,k+1) - gyz(i,j,k-1))/(2.*cctk_delta_space(3))

    ! gzz 
    gzzdx(i,j,k) = (gzz(i+1,j,k) - gzz(i-1,j,k))/(2.*cctk_delta_space(1))
    gzzdy(i,j,k) = (gzz(i,j+1,k) - gzz(i,j-1,k))/(2.*cctk_delta_space(2))
    gzzdz(i,j,k) = (gzz(i,j,k+1) - gzz(i,j,k-1))/(2.*cctk_delta_space(3))

end subroutine central_finite