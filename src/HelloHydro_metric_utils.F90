#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Collecting of subroutines for calculate metric quanties that are used by
! a variety of routines
! Should NOT be scheduled by Cactus, and only accessed by a scheduled routine 
module metric_utils
    implicit none 
    contains 
    subroutine get_covshift(CCTK_ARGUMENTS,i,j,k,betacox,betacoy,betacoz)
        implicit none
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_FUNCTIONS
        DECLARE_CCTK_PARAMETERS
        integer, intent(in) :: i,j,k
        real*8, intent(out) :: betacox, betacoy, betacoz

        !print*, "In covshift"
        betacox = gxx(i,j,k)*betax(i,j,k) + gxy(i,j,k)*betay(i,j,k) + gxz(i,j,k)*betaz(i,j,k)
        betacoy = gxy(i,j,k)*betax(i,j,k) + gyy(i,j,k)*betay(i,j,k) + gyz(i,j,k)*betaz(i,j,k)
        betacoz = gxz(i,j,k)*betax(i,j,k) + gyz(i,j,k)*betay(i,j,k) + gzz(i,j,k)*betaz(i,j,k)
 
    end subroutine get_covshift

    subroutine get_dtcovshift(CCTK_ARGUMENTS,i,j,k,dtbetacox,dtbetacoy,dtbetacoz)
        implicit none
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_FUNCTIONS
        DECLARE_CCTK_PARAMETERS
        integer, intent(in) :: i,j,k
        real*8, intent(out) :: dtbetacox, dtbetacoy, dtbetacoz

        dtbetacox = gxx(i,j,k)*dtbetax(i,j,k) + gxy(i,j,k)*dtbetay(i,j,k) + gxz(i,j,k)*dtbetaz(i,j,k)
        dtbetacoy = gxy(i,j,k)*dtbetax(i,j,k) + gyy(i,j,k)*dtbetay(i,j,k) + gyz(i,j,k)*dtbetaz(i,j,k)
        dtbetacoz = gxz(i,j,k)*dtbetax(i,j,k) + gyz(i,j,k)*dtbetay(i,j,k) + gzz(i,j,k)*dtbetaz(i,j,k)
 
    end subroutine get_dtcovshift

    subroutine HelloHydro_metric_global(CCTK_ARGUMENTS)
        use einsteintk_utils
        DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_FUNCTIONS
        DECLARE_CCTK_PARAMETERS
        integer :: i,j,k,ierr
        integer :: il,jl,kl
    
        ! Probably really slow to recreate these large grids
        ! everytime we want to pass to phantom
        ! These should be initialized and then left
        ! Better yet, use the arrays from einsteintk_utils in phantom
        ! i.e call phantom init which sets up the arrays and then set it with the global processor
        ! CCTK_REAL, allocatable  :: gcovgrid(:,:,:,:,:)
        ! CCTK_REAL, allocatable  :: gcongrid(:,:,:,:,:)
        ! CCTK_REAL, allocatable  :: sqrtggrid(:,:,:)
        ! CCTK_REAL, allocatable  :: tmunugrid(:,:,:,:,:)
        ! CCTK_REAL, allocatable  :: posgrid(:,:,:) ! The position of the grid (centre??) for a given x,y,z val
    
        ! !If we are on the main processor
        ! if (CCTK_MyProc(cctkGH)==0) then
        !     call CCTK_INFO('Setting up global arrays (for testing only!!)')
        !     ! Again don't need this but just for testing 
        !     allocate(gcovgrid(4,4,cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)))
        !     allocate(gcongrid(4,4,cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)))
        !     allocate(sqrtggrid(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)))
        !     allocate(posgrid(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)))
        !     allocate(tmunugrid(4,4,cctk_gsh(1),cctk_gsh(2),cctk_gsh(3)))
        !     call CCTK_INFO('Done setting up global arrays!')
        ! endif 
    
        ! ! All other processors wait here until we've allocated the arrays 
        ! call CCTK_Barrier(ierr,cctkGH) 
    
        ! Get the position (index) of the local grid in the global grid
        il = cctk_lbnd(1)! + 1
        jl = cctk_lbnd(2)! + 1  ! indices output from cctk_lbnd start at 0, need to +1
        kl = cctk_lbnd(3)! + 1

        print*, "i lower: ", il
        print*, "j lower: ", jl
        print*, "k lower: ", kl
    
        ! Copy local values to the global grid 
        do k=1, cctk_lsh(3)
            do j=1, cctk_lsh(2)
                do i=1, cctk_lsh(1)
    
                    ! Check that this is correct
                    ! Set the covariant grid
                    ! Gcovgrid runs from 0-3 as this is the convention used in phantom for metrics 
                    gcovgrid(0,0,i + il, j + jl, k + kl) = gmunutt(i,j,k)
                    !if (i==1 .and. j==1 .and. k==1) print*, "gmunutt inner loop: ", gcovgrid(i,j,k)
                    gcovgrid(0,1,i + il, j + jl, k + kl) = gmunutx(i,j,k)
                    gcovgrid(0,2,i + il, j + jl, k + kl) = gmunuty(i,j,k)
                    gcovgrid(0,3,i + il, j + jl, k + kl) = gmunutz(i,j,k)
                    gcovgrid(1,0,i + il, j + jl, k + kl) = gmunutx(i,j,k)
                    gcovgrid(1,1,i + il, j + jl, k + kl) = gmunuxx(i,j,k)
                    gcovgrid(1,2,i + il, j + jl, k + kl) = gmunuxy(i,j,k)
                    gcovgrid(1,3,i + il, j + jl, k + kl) = gmunuxz(i,j,k)
                    gcovgrid(2,0,i + il, j + jl, k + kl) = gmunuty(i,j,k)
                    gcovgrid(2,1,i + il, j + jl, k + kl) = gmunuxy(i,j,k)
                    gcovgrid(2,2,i + il, j + jl, k + kl) = gmunuyy(i,j,k)
                    gcovgrid(2,3,i + il, j + jl, k + kl) = gmunuyz(i,j,k)
                    gcovgrid(3,0,i + il, j + jl, k + kl) = gmunutz(i,j,k)
                    gcovgrid(3,1,i + il, j + jl, k + kl) = gmunuxz(i,j,k)
                    gcovgrid(3,2,i + il, j + jl, k + kl) = gmunuyz(i,j,k)
                    gcovgrid(3,3,i + il, j + jl, k + kl) = gmunuzz(i,j,k)
                    
                    ! Set the contravariant grid 
                    gcongrid(0,0,i + il, j + jl, k + kl) = gmunucontt(i,j,k)
                    gcongrid(0,1,i + il, j + jl, k + kl) = gmunucontx(i,j,k)
                    gcongrid(0,2,i + il, j + jl, k + kl) = gmunuconty(i,j,k)
                    gcongrid(0,3,i + il, j + jl, k + kl) = gmunucontz(i,j,k)
                    gcongrid(1,0,i + il, j + jl, k + kl) = gmunucontx(i,j,k)
                    gcongrid(1,1,i + il, j + jl, k + kl) = gmunuconxx(i,j,k)
                    gcongrid(1,2,i + il, j + jl, k + kl) = gmunuconxy(i,j,k)
                    gcongrid(1,3,i + il, j + jl, k + kl) = gmunuconxz(i,j,k)
                    gcongrid(2,0,i + il, j + jl, k + kl) = gmunuconty(i,j,k)
                    gcongrid(2,1,i + il, j + jl, k + kl) = gmunuconxy(i,j,k)
                    gcongrid(2,2,i + il, j + jl, k + kl) = gmunuconyy(i,j,k)
                    gcongrid(2,3,i + il, j + jl, k + kl) = gmunuconyz(i,j,k)
                    gcongrid(3,0,i + il, j + jl, k + kl) = gmunucontz(i,j,k)
                    gcongrid(3,1,i + il, j + jl, k + kl) = gmunuconxz(i,j,k)
                    gcongrid(3,2,i + il, j + jl, k + kl) = gmunuconyz(i,j,k)
                    gcongrid(3,3,i + il, j + jl, k + kl) = gmunuconzz(i,j,k)
    
                    ! Set neg sqrtg 
                    sqrtggrid(i + il, j + jl, k + kl) = sqrtg(i,j,k)
    
                enddo 
            enddo 
        enddo 
        
        ! This is to test that we are actually initialising some values
        print*, "gcovgrid: ", gcovgrid(1,1,1, 1, 1)
        print*, "gmunutt: ", gmunutt(1,1,1)
        print*, "Local grid sizes are : ", cctk_lsh(1), cctk_lsh(2), cctk_lsh(3)
        print*, "Global grid size is: ", cctk_gsh(1), cctk_gsh(2), cctk_gsh(3)
        !stop 
    end subroutine HelloHydro_metric_global
end module metric_utils