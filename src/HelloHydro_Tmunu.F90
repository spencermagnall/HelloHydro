#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating the stress energy tensor and adding its components

! Eventually this should be similar to GRHydro_Tmunu and calculate Tmunu properly but for now we just use dust solution

! I assume that everything is MCRF and that the 4 velocity is just u0 = 1, but this could be incorrect?? 
subroutine HelloHydro_Tmunu(CCTK_ARGUMENTS)
    use einsteintk_utils
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer :: i,j,k
    call CCTK_INFO("Calculating components of stress energy tensor from primative variables")
    
    !print*, "Value of Stress energy tensor before is: ", eTtt(6,6,6)
    !print*, "tmunugrid: ", tmunugrid(1,1,6,6,6)
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                 

                    ! I really hate the way ET splits up its variables 
                ! This is super confusing using zero based indexing in
                ! Phantom but not in ET
                ! I Need to think of a better way to do this 
                ! As I should not directly be setting the Tmunu
                ! Scalar part of Tmunu               
                eTtt(i,j,k) = tmunugrid(0,0,i,j,k)
               
                ! Vector part of Tmunu
                eTtx(i,j,k) = tmunugrid(0,1,i,j,k) 
                eTty(i,j,k) = tmunugrid(0,2,i,j,k)
                eTtz(i,j,k) = tmunugrid(0,3,i,j,k) 
                ! Tensor part ov Tmunu
                eTxx(i,j,k) = tmunugrid(1,1,i,j,k)
                eTxy(i,j,k) = tmunugrid(1,2,i,j,k)
                eTxz(i,j,k) = tmunugrid(1,3,i,j,k)
                eTyy(i,j,k) = tmunugrid(2,2,i,j,k)
                eTyz(i,j,k) = tmunugrid(2,3,i,j,k)
                eTzz(i,j,k) = tmunugrid(3,3,i,j,k)

            enddo 
        enddo 
    enddo 
    !print*, "Value of Stress energy tensor is: ", eTtt(1,1,1)
    print*,"Stress energy tensor: ", eTtt(6,6,6), eTxx(6,6,6), eTyy(6,6,6), eTzz(6,6,6)
    print*, "values of metric are: ", gxx(1,1,1), gyy(1,1,1), gzz(1,1,1)
    
    !stop 
end subroutine HelloHydro_Tmunu
