#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

! Collecting of subroutines for calculate metric quanties that are used by
! a variety of routines
! Should NOT be scheduled by Cactus, and only accessed by a scheduled routine 
subroutine get_covshift(CCTK_ARGUMENTS,i,j,k,betacox,betacoy,betacoz)
    implicit none
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    integer, intent(in) :: i,j,k
    real*8, intent(out) :: betacox, betacoy, betacoz

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