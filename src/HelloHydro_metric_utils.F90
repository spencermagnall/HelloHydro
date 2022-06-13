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

                    ! Set metric derivs
                    ! gttdx
                    metricderivsgrid(0,0,1,i + il, j + jl, k + kl) = gttdx(i,j,k)
                    ! gttdy 
                    metricderivsgrid(0,0,2,i + il, j + jl, k + kl) = gttdy(i,j,k)
                    ! gttdz
                    metricderivsgrid(0,0,3,i + il, j + jl, k + kl) = gttdz(i,j,k)
                    
                    ! gtxdx
                    metricderivsgrid(0,1,1,i + il, j + jl, k + kl) = gtxdx(i,j,k)
                    ! gtxdy 
                    metricderivsgrid(0,1,2,i + il, j + jl, k + kl) = gtxdy(i,j,k)
                    ! gtxdz
                    metricderivsgrid(0,1,3,i + il, j + jl, k + kl) = gtxdz(i,j,k)

                    ! Copy to the symmetrical comp gxt
                    metricderivsgrid(1,0,:,i + il, j + jl, k + kl) = metricderivsgrid(0,1,:,i + il, j + jl, k + kl) 

                    ! gtydx
                    metricderivsgrid(0,2,1,i + il, j + jl, k + kl) = gtydx(i,j,k)
                    ! gtydy 
                    metricderivsgrid(0,2,2,i + il, j + jl, k + kl) = gtydy(i,j,k)
                    ! gtydz
                    metricderivsgrid(0,2,3,i + il, j + jl, k + kl) = gtydz(i,j,k)
                    
                    ! Copy to the symmetrical comp gyt
                    metricderivsgrid(2,0,:,i + il, j + jl, k + kl) = metricderivsgrid(0,2,:,i + il, j + jl, k + kl) 

                    ! gtzdx
                    metricderivsgrid(0,3,1,i + il, j + jl, k + kl) = gtzdx(i,j,k)
                    ! gtzdy 
                    metricderivsgrid(0,3,2,i + il, j + jl, k + kl) = gtzdy(i,j,k)
                    ! gtzdz
                    metricderivsgrid(0,3,3,i + il, j + jl, k + kl) = gtzdz(i,j,k)

                    ! Copy to the symmetrical comp gzt
                    metricderivsgrid(3,0,:,i + il, j + jl, k + kl) = metricderivsgrid(0,3,:,i + il, j + jl, k + kl) 

                    ! gxxdx
                    metricderivsgrid(1,1,1,i + il, j + jl, k + kl) = gxxdx(i,j,k)
                    ! gxxdy 
                    metricderivsgrid(1,1,2,i + il, j + jl, k + kl) = gxxdy(i,j,k)
                    ! gxxdz
                    metricderivsgrid(1,1,3,i + il, j + jl, k + kl) = gxxdz(i,j,k)

                    ! gxydx
                    metricderivsgrid(1,2,1,i + il, j + jl, k + kl) = gxydx(i,j,k)
                    ! gxydy 
                    metricderivsgrid(1,2,2,i + il, j + jl, k + kl) = gxydy(i,j,k)
                    ! gxydz
                    metricderivsgrid(1,2,3,i + il, j + jl, k + kl) = gxydz(i,j,k)

                    ! Copy to the symmetrical comp gyx
                    metricderivsgrid(2,1,:,i + il, j + jl, k + kl) = metricderivsgrid(1,2,:,i + il, j + jl, k + kl)! gtzdx
                    
                    ! gxzdx
                    metricderivsgrid(1,3,1,i + il, j + jl, k + kl) = gxzdx(i,j,k)
                    ! gxzdy 
                    metricderivsgrid(1,3,2,i + il, j + jl, k + kl) = gxzdy(i,j,k)
                    ! gxzdz
                    metricderivsgrid(1,3,3,i + il, j + jl, k + kl) = gxzdz(i,j,k)

                    ! Copy to the symmetrical comp gzx
                    metricderivsgrid(3,1,:,i + il, j + jl, k + kl) = metricderivsgrid(3,1,:,i + il, j + jl, k + kl) 

                    ! gyydx
                    metricderivsgrid(2,2,1,i + il, j + jl, k + kl) = gyydx(i,j,k)
                    ! gyydy 
                    metricderivsgrid(2,2,2,i + il, j + jl, k + kl) = gyydy(i,j,k)
                    ! gyydz
                    metricderivsgrid(2,2,3,i + il, j + jl, k + kl) = gyydz(i,j,k)

                    ! gyzdx
                    metricderivsgrid(2,3,1,i + il, j + jl, k + kl) = gyzdx(i,j,k)
                    ! gyzdy 
                    metricderivsgrid(2,3,2,i + il, j + jl, k + kl) = gyzdy(i,j,k)
                    ! gyzdz
                    metricderivsgrid(2,3,3,i + il, j + jl, k + kl) = gyzdz(i,j,k)

                    ! Copy to the symmetrical comp gzt
                    metricderivsgrid(3,2,:,i + il, j + jl, k + kl) = metricderivsgrid(2,3,:,i + il, j + jl, k + kl) 
                    
                    ! gzzdx
                    metricderivsgrid(3,3,1,i + il, j + jl, k + kl) = gzzdx(i,j,k)
                    ! gzzdy 
                    metricderivsgrid(3,3,2,i + il, j + jl, k + kl) = gzzdy(i,j,k)
                    ! gzzdz
                    metricderivsgrid(3,3,3,i + il, j + jl, k + kl) = gzzdz(i,j,k)


    
                enddo 
            enddo 
        enddo 
        
        ! This is to test that we are actually initialising some values
        print*, "gcovgrid: ", gcovgrid(:,:,1, 1, 1)
        print*, "gmunutt: ", gmunutt(1,1,1)
        print*, "Local grid sizes are : ", cctk_lsh(1), cctk_lsh(2), cctk_lsh(3)
        print*, "Global grid size is: ", cctk_gsh(1), cctk_gsh(2), cctk_gsh(3)
        print*, "sqrtggrid: ", sqrtggrid(1,1,1)
        !stop

    end subroutine HelloHydro_metric_global



end module metric_utils