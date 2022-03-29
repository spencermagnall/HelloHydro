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
    integer :: i,j,k
    call CCTK_INFO("Calculating components of stress energy tensor from primative variables")
    
    print*, "Value of Stress energy tensor before is: ", eTtt(1,1,1)
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                 

                    eTtt(i,j,k) = eTtt(i,j,k) + rho(i,j,k)
                    eTxx(i,j,k) = eTxx(i,j,k) + press(i,j,k)
                    eTyy(i,j,k) = eTyy(i,j,k) + press(i,j,k)
                    eTzz(i,j,k) = eTzz(i,j,k) + press(i,j,k)

            enddo 
        enddo 
    enddo 
    !print*, "Value of Stress energy tensor is: ", eTtt(1,1,1)
    print*,"Stress energy tensor: ", eTtt(1,1,1), eTtx(1,1,1), eTty(1,1,    1), eTtz(1,1,1), eTxx(1,1,1), eTyy(1,1,1), eTzz(1,1,1), &
         eTxy(1,1,1), eTxz(1,1,1), eTyz(1,1,1)
    print*, "values of metric are: ", gxx(1,1,1), gyy(1,1,1), gzz(1,1,1)
    
    !stop 
end subroutine HelloHydro_Tmunu
