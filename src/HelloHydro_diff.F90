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
    
    print*, "Spatial derivs: ", gxxdx(5,5,5),gxxdy(5,5,5),gxxdz(1,1,1)
    
end subroutine HelloHydro_diff

subroutine central_finite(CCTK_ARGUMENTS,i,j,k)
    implicit none 
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: i,j,k
    integer             :: iupper, jupper, kupper
    integer             :: ilower, jlower, klower 

    
    ! Get first derivatives of the spatial metric using a centered finite difference 
    ! Refactor later into a new function as there is a lot of copy and paste here 
    ! Might also be a better idea to use a general difff class to alow the code to be more modular 
    
    ! Periodic boundary conditions are hanndled here
    ! Again this should be in its own subroutine 
    
    call get_lower_upper_index(CCTK_ARGUMENTS,i,1,ilower,iupper)
    call get_lower_upper_index(CCTK_ARGUMENTS,j,2,jlower,jupper)
    call get_lower_upper_index(CCTK_ARGUMENTS,k,3,klower,kupper)
   
    ! gxx 
    gxxdx(i,j,k) = (gxx(iupper,j,k) - gxx(ilower,j,k))/(2.*cctk_delta_space(1))
    !print*, "gxx values: ", gxx(i+1,j,k), gxx(i-1,j,k), gxxdx(i,j,k)
    !stop
    gxxdy(i,j,k) = (gxx(i,jupper,k) - gxx(i,jlower,k))/(2.*cctk_delta_space(2))
    gxxdz(i,j,k) = (gxx(i,j,kupper) - gxx(i,j,klower))/(2.*cctk_delta_space(3))
    ! gxy 
    gxydx(i,j,k) = (gxy(iupper,j,k) - gxy(ilower,j,k))/(2.*cctk_delta_space(1))
    gxydy(i,j,k) = (gxy(i,jupper,k) - gxy(i,jlower,k))/(2.*cctk_delta_space(2))
    gxydz(i,j,k) = (gxy(i,j,kupper) - gxy(i,j,klower))/(2.*cctk_delta_space(3))
    ! gxz
    gxzdx(i,j,k) = (gxz(iupper,j,k) - gxz(ilower,j,k))/(2.*cctk_delta_space(1))
    gxzdy(i,j,k) = (gxz(i,jupper,k) - gxz(i,jlower,k))/(2.*cctk_delta_space(2))
    gxzdz(i,j,k) = (gxz(i,j,kupper) - gxz(i,j,klower))/(2.*cctk_delta_space(3))
    ! gyy 
    gyydx(i,j,k) = (gyy(iupper,j,k) - gyy(ilower,j,k))/(2.*cctk_delta_space(1))
    gyydy(i,j,k) = (gyy(i,jupper,k) - gyy(i,jlower,k))/(2.*cctk_delta_space(2))
    gyydz(i,j,k) = (gyy(i,j,kupper) - gyy(i,j,klower))/(2.*cctk_delta_space(3))

    ! gyz 
    gyzdx(i,j,k) = (gyz(iupper,j,k) - gyz(ilower,j,k))/(2.*cctk_delta_space(1))
    gyzdy(i,j,k) = (gyz(i,jupper,k) - gyz(i,jlower,k))/(2.*cctk_delta_space(2))
    gyzdz(i,j,k) = (gyz(i,j,kupper) - gyz(i,j,klower))/(2.*cctk_delta_space(3))

    ! gzz 
    gzzdx(i,j,k) = (gzz(iupper,j,k) - gzz(ilower,j,k))/(2.*cctk_delta_space(1))
    gzzdy(i,j,k) = (gzz(i,jupper,k) - gzz(i,jlower,k))/(2.*cctk_delta_space(2))
    gzzdz(i,j,k) = (gzz(i,j,kupper) - gzz(i,j,klower))/(2.*cctk_delta_space(3))

end subroutine central_finite

subroutine get_lower_upper_index(CCTK_ARGUMENTS,index,dir,indexlower,indexupper)
    implicit none 
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: index, dir ! Index to get upper and lower bound for and direction i.e x
    integer, intent(out) :: indexlower, indexupper
    
    if (index == 1) then  
        indexlower = cctk_lsh(dir)
    else 
        indexlower = index - 1 
    endif
    if (index == cctk_lsh(dir)) then  
        indexupper = 1
    else 
        indexupper = index + 1
    endif  
end subroutine get_lower_upper_index