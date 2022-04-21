#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Subroutine for calculating metric quantites 
subroutine HelloHydro_metric(CCTK_ARGUMENTS)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer :: i,j,k
    call CCTK_INFO("Calculating metric quantities for the grid ")
    
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                 
                ! Reconstruct 4 metric from the lapse and shift and the spatial metric 
                call get_4metric(CCTK_ARGUMENTS,i,j,k)
                call get_sqrtg(CCTK_ARGUMENTS,i,j,k)
                call get_con(CCTK_ARGUMENTS,i,j,k)
                ! get metric derivatives through finite difference
                ! Which finite difference should we use?
                ! How does this work with mpi do we use ghost points?
            enddo 
        enddo 
    enddo 
    print*, "gmunutt, gmunuxx, gmunuyy, gmunuzz ", gmunutt(1,1,1), gmunuxx(1,1,1), gmunuyy(1,1,1),gmunuzz(1,1,1)
    print*, "Non-diagonial terms: ", gmunutx(1,1,1),gmunuty(1,1,1),gmunutz(1,1,1),gmunuxy(1,1,1),gmunuxz(1,1,1),gmunuyz(1,1,1)
    print*, "gmunucontt, gmunuconxx, gmunuconyy, gmunuconzz ", gmunucontt(1,1,1), gmunuconxx(1,1,1), gmunuconyy(1,1,1),gmunuconzz(1,1,1)
    print*, "Non-diagonial terms: ", gmunucontx(1,1,1),gmunuconty(1,1,1),gmunucontz(1,1,1),gmunuconxy(1,1,1),gmunuconxz(1,1,1),gmunuconyz(1,1,1)
    print*, "Sqrtg: ", sqrtg(1,1,1)
    !stop 
end subroutine HelloHydro_metric

subroutine get_4metric(CCTK_ARGUMENTS,i,j,k)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: i,j,k
    real*8 ::  beta2

    ! Get the dot product of the shift vector
    beta2 = betax(i,j,k)*betax(i,j,k)  + betay(i,j,k)*betay(i,j,k) + betaz(i,j,k)*betaz(i,j,k)
    !print*, "Beta squared: ", beta2
    !print*, "gmunutt: ", gmunutt
    !print*, "gmunutt: ", gmunutt(i,j,k)
    !print*, "rhs: ", -(alp(i,j,k)*alp(i,j,k) - beta2)
    !stop
    gmunutt(i,j,k) = -(alp(i,j,k)*alp(i,j,k) - beta2)
    !stop 
    gmunutx(i,j,k) = 2.*betax(i,j,k)
    gmunuty(i,j,k) = 2.*betay(i,j,k)
    gmunutz(i,j,k) = 2.*betaz(i,j,k)
    ! Note the use of symmetry of the metric tensor 
    gmunuxx(i,j,k) = gxx(i,j,k)
    gmunuxy(i,j,k) = gxy(i,j,k)
    gmunuxz(i,j,k) = gxz(i,j,k)
    gmunuyy(i,j,k) = gyy(i,j,k)
    gmunuyz(i,j,k) = gyz(i,j,k)
    gmunuzz(i,j,k) = gzz(i,j,k)

end subroutine get_4metric

subroutine get_con(CCTK_ARGUMENTS,i,j,k)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: i,j,k
    real*8 :: gmunu(4,4), gmunucon(4,4)
    real*8 :: det, det1 
    real*8 :: a11,a12,a13,a14
    real*8 :: a21,a22,a23,a24
    real*8 :: a31,a32,a33,a34
    real*8 :: a41,a42,a43,a44
    
    ! This is probably a bit slow but makes the code easier
    ! Generate local copy of gmunnu
    call get_local_gmunu(CCTK_ARGUMENTS,i,j,k,gmunu)

    a11 = gmunu(1,1)
    a21 = gmunu(2,1)
    a31 = gmunu(3,1)
    a41 = gmunu(4,1)
    a12 = gmunu(1,2)
    a22 = gmunu(2,2)
    a32 = gmunu(3,2)
    a42 = gmunu(4,2)
    a13 = gmunu(1,3)
    a23 = gmunu(2,3)
    a33 = gmunu(3,3)
    a43 = gmunu(4,3)
    a14 = gmunu(1,4)
    a24 = gmunu(2,4)
    a34 = gmunu(3,4)
    a44 = gmunu(4,4)


        
    call get_det(gmunu,det)

    det1 = 1./det

    gmunucon(1,1) = -(a24*a33*a42) + a23*a34*a42 + a24*a32*a43 - a22*a34*a43 - a23*a32*   a44 + a22*a33*a44
    gmunucon(2,1) =    a14*a33*a42 - a13*a34*a42 - a14*a32*a43 + a12*a34*a43 + a13*a32*   a44 - a12*a33*a44
    gmunucon(3,1) = -(a14*a23*a42) + a13*a24*a42 + a14*a22*a43 - a12*a24*a43 - a13*a22*   a44 + a12*a23*a44
    gmunucon(4,1) =    a14*a23*a32 - a13*a24*a32 - a14*a22*a33 + a12*a24*a33 + a13*a22*   a34 - a12*a23*a34
    gmunucon(1,2) =    a24*a33*a41 - a23*a34*a41 - a24*a31*a43 + a21*a34*a43 + a23*a31*   a44 - a21*a33*a44
    gmunucon(2,2) = -(a14*a33*a41) + a13*a34*a41 + a14*a31*a43 - a11*a34*a43 - a13*a31*   a44 + a11*a33*a44
    gmunucon(3,2) =    a14*a23*a41 - a13*a24*a41 - a14*a21*a43 + a11*a24*a43 + a13*a21*   a44 - a11*a23*a44
    gmunucon(4,2) = -(a14*a23*a31) + a13*a24*a31 + a14*a21*a33 - a11*a24*a33 - a13*a21*   a34 + a11*a23*a34
    gmunucon(1,3) = -(a24*a32*a41) + a22*a34*a41 + a24*a31*a42 - a21*a34*a42 - a22*a31*   a44 + a21*a32*a44
    gmunucon(2,3) =    a14*a32*a41 - a12*a34*a41 - a14*a31*a42 + a11*a34*a42 + a12*a31*   a44 - a11*a32*a44
    gmunucon(3,3) = -(a14*a22*a41) + a12*a24*a41 + a14*a21*a42 - a11*a24*a42 - a12*a21*   a44 + a11*a22*a44
    gmunucon(4,3) =    a14*a22*a31 - a12*a24*a31 - a14*a21*a32 + a11*a24*a32 + a12*a21*   a34 - a11*a22*a34
    gmunucon(1,4) =    a23*a32*a41 - a22*a33*a41 - a23*a31*a42 + a21*a33*a42 + a22*a31*   a43 - a21*a32*a43
    gmunucon(2,4) = -(a13*a32*a41) + a12*a33*a41 + a13*a31*a42 - a11*a33*a42 - a12*a31*   a43 + a11*a32*a43
    gmunucon(3,4) =    a13*a22*a41 - a12*a23*a41 - a13*a21*a42 + a11*a23*a42 + a12*a21*   a43 - a11*a22*a43
    gmunucon(4,4) = -(a13*a22*a31) + a12*a23*a31 + a13*a21*a32 - a11*a23*a32 - a12*a21*   a33 + a11*a22*a33
    gmunucon = gmunucon*det1
    
    ! Set the cactus variables 
    gmunucontt(i,j,k) = gmunucon(1,1)
    gmunucontx(i,j,k) = gmunucon(1,2)
    gmunuconty(i,j,k) = gmunucon(1,3)
    gmunucontz(i,j,k) = gmunucon(1,4)
    gmunuconxx(i,j,k) = gmunucon(2,2)
    gmunuconxy(i,j,k) = gmunucon(2,3)
    gmunuconxz(i,j,k) = gmunucon(2,4)
    gmunuconyy(i,j,k) = gmunucon(3,3)
    gmunuconyz(i,j,k) = gmunucon(3,4)
    gmunuconzz(i,j,k) = gmunucon(4,4)

end subroutine get_con

subroutine get_sqrtg(CCTK_ARGUMENTS,i,j,k)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: i,j,k
    real*8 :: gmunu(4,4)
    !real, intent(out) :: sqrtg
    real*8 :: det
   
    ! TODO What about when g is negative as for the flrw metric? how do we handle that case 
    ! We want sqrt of -det not sqrt of det
    !Get local copy of gmunu
    call get_local_gmunu(CCTK_ARGUMENTS,i,j,k,gmunu)
    ! Get det
    call get_det(gmunu,det)
    !print*, "det: ",det
    sqrtg(i,j,k) = sqrt(-det)
end subroutine get_sqrtg

subroutine get_det(gmunu,det)
    real*8, intent(in)  :: gmunu(4,4)
    real*8, intent(out) :: det
    real*8 :: a11,a12,a13,a14
    real*8 :: a21,a22,a23,a24
    real*8 :: a31,a32,a33,a34
    real*8 :: a41,a42,a43,a44

    a11 = gmunu(1,1)
    a21 = gmunu(2,1)
    a31 = gmunu(3,1)
    a41 = gmunu(4,1)
    a12 = gmunu(1,2)
    a22 = gmunu(2,2)
    a32 = gmunu(3,2)
    a42 = gmunu(4,2)
    a13 = gmunu(1,3)
    a23 = gmunu(2,3)
    a33 = gmunu(3,3)
    a43 = gmunu(4,3)
    a14 = gmunu(1,4)
    a24 = gmunu(2,4)
    a34 = gmunu(3,4)
    a44 = gmunu(4,4)

    det = a14*a23*a32*a41 - a13*a24*a32*a41 - a14*a22*a33*a41 + a12*a24*a33*a41 +  &
        a13*a22*a34*a41 - a12*a23*a34*a41 - a14*a23*a31*a42 + a13*a24*a31*a42 +  &
        a14*a21*a33*a42 - a11*a24*a33*a42 - a13*a21*a34*a42 + a11*a23*a34*a42 +  &
        a14*a22*a31*a43 - a12*a24*a31*a43 - a14*a21*a32*a43 + a11*a24*a32*a43 +  &
        a12*a21*a34*a43 - a11*a22*a34*a43 - a13*a22*a31*a44 + a12*a23*a31*a44 +  &
        a13*a21*a32*a44 - a11*a23*a32*a44 - a12*a21*a33*a44 + a11*a22*a33*a44

end subroutine get_det

subroutine get_local_gmunu(CCTK_ARGUMENTS,i,j,k,gmunu)
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
 
    integer, intent(in) :: i,j,k
    real*8, intent(out)   :: gmunu(4,4)

    gmunu(1,1) = gmunutt(i,j,k)
    gmunu(1,2) = gmunutx(i,j,k)
    gmunu(1,3) = gmunuty(i,j,k)
    gmunu(1,4) = gmunutz(i,j,k)
    gmunu(2,1) = gmunutx(i,j,k)
    gmunu(2,2) = gmunuxx(i,j,k)
    gmunu(2,3) = gmunuxy(i,j,k)
    gmunu(2,4) = gmunuxz(i,j,k)
    gmunu(3,1) = gmunuty(i,j,k)
    gmunu(3,2) = gmunuxy(i,j,k)
    gmunu(3,3) = gmunuyy(i,j,k)
    gmunu(3,4) = gmunuyz(i,j,k)
    gmunu(4,1) = gmunutz(i,j,k)
    gmunu(4,2) = gmunuxz(i,j,k)
    gmunu(4,3) = gmunuyz(i,j,k)
    gmunu(4,4) = gmunuzz(i,j,k)

 
end subroutine get_local_gmunu
