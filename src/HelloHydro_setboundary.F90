#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating setting the stress-energy tensor for ghost cells which
! Lie outside of the physical domain 
! Only valid for periodic boundary conditions

subroutine HelloHydro_setboundary_Tmunu(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer :: i,j,k
    integer :: nghostx, nghosty, nghostz
    integer :: indexghostx, indexghosty, indexghostz ! The index of the ghost zones 
    integer :: indexcopyx, indexcopyy, indexcopyz    ! The index which the stress energy tensor is copied from
    call CCTK_INFO("Setting the stress energy for the boundary points")

    ! Work out the number of ghost points in each direction
    nghostx = cctk_nghostzones(1)
    nghosty = cctk_nghostzones(2) 
    nghostz = cctk_nghostzones(3)

    ! Copy the stress energy tensor for the lower ghost zones 
    ! We want to go from lower grid bound to nghostx 

    do k=1, nghostz
        do j=1, nghosty
            do i=1, nghostx

                ! Copy the boundary conditions from the other side (assuming periodicity)
                ! Note that we DO NOT want the ghost zones on the values from the ghost
                ! zones on the other side just the interior quantities
                
                indexcopyx = cctk_gsh(1) - nghostx - i 
                indexcopyy = cctk_gsh(2) - nghosty - j 
                indexcopyz = cctk_gsh(3) - nghostz - k

                ! Scalar part of Tmunu               
                eTtt(i,j,k) = eTtt(indexcopyx, indexcopyy, indexcopyz)
               
                ! Vector part of Tmunu
                eTtx(i,j,k) =  eTtx(indexcopyx, indexcopyy, indexcopyz) 
                eTty(i,j,k) =  eTty(indexcopyx, indexcopyy, indexcopyz)
                eTtz(i,j,k) =  eTtz(indexcopyx, indexcopyy, indexcopyz) 
                ! Tensor part ov Tmunu
                eTxx(i,j,k) =  eTxx(indexcopyx, indexcopyy, indexcopyz)
                eTxy(i,j,k) =  eTxy(indexcopyx, indexcopyy, indexcopyz)
                eTxz(i,j,k) =  eTxz(indexcopyx, indexcopyy, indexcopyz)
                eTyy(i,j,k) =  eTyy(indexcopyx, indexcopyy, indexcopyz)
                eTyz(i,j,k) =  eTyz(indexcopyx, indexcopyy, indexcopyz)
                eTzz(i,j,k) =  eTzz(indexcopyx, indexcopyy, indexcopyz)

            enddo 
        enddo 
    enddo 


    ! Copy the stress energy tensor for the upper ghost zones
    ! We want to go from upperbound - nghostx to upperbound

    do k=1, nghostz
        do j=1, nghosty
            do i=1, nghostx
               
                indexcopyx = nghostx + i 
                indexcopyy = nghosty + j 
                indexcopyz = nghostz + k 

                indexghostx = cctk_gsh(1)-nghostx-i
                indexghosty = cctk_gsh(2)-nghosty-j 
                indexghostz = cctk_gsh(3)-nghostz-k

                ! Copy the boundary conditions from the other side
                ! Scalar part of Tmunu               
                eTtt(indexghostx, indexghosty, indexghostz) = eTtt(indexcopyx, indexcopyy, indexcopyz)
               
                ! Vector part of Tmunu
                eTtx(indexghostx, indexghosty, indexghostz)   = eTtx(indexcopyx, indexcopyy, indexcopyz)
                eTty(indexghostx, indexghosty, indexghostz)   = eTty(indexcopyx, indexcopyy, indexcopyz) 
                eTtz(indexghostx, indexghosty, indexghostz)   = eTtz(indexcopyx, indexcopyy, indexcopyz) 
                ! Tensor part ov Tmunu
                eTxx(indexghostx, indexghosty, indexghostz)   = eTxx(indexcopyx, indexcopyy, indexcopyz)
                eTxy(indexghostx,indexghosty,indexghostz) =  eTxy(indexcopyx, indexcopyy, indexcopyz)
                eTxz(indexghostx,indexghosty,indexghostz) =  eTxz(indexcopyx, indexcopyy, indexcopyz)
                eTyy(indexghostx,indexghosty,indexghostz) =  eTyy(indexcopyx, indexcopyy, indexcopyz)
                eTyz(indexghostx,indexghosty,indexghostz) =  eTyz(indexcopyx, indexcopyy, indexcopyz)
                eTzz(indexghostx,indexghosty,indexghostz) =  eTzz(indexcopyx, indexcopyy, indexcopyz) 
            enddo 
        enddo 
    enddo  
    print*, "eTtt bounds: ", eTtt(1:3,1:3,1:3)
    stop 
end subroutine HelloHydro_setboundary_Tmunu
