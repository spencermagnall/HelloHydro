#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine HelloHydro_init(CCTK_ARGUMENTS)
    use io,              only:id,master,nprocs,set_io_unit_numbers,die
    use initial,         only:initialise,finalise,startrun,endrun
    implicit none
    character(len=120) :: infile,logfile,evfile,dumpfile,path
    integer :: i,j,k,pathstringlength
    real :: poly_k, poly_gamma
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    
    poly_k = 1.e-20
    poly_gamma = 2.0
    ! For now we just hardcode the infile, to see if startrun actually works!
    ! I'm not sure what the best way to actually do this is? 
    ! Do we store the phantom.in file in par and have it read from there?
    !infile = "/Users/spencer/phantomET/phantom/test/flrw.in"
    !infile = trim(infile)//'.in'
    !print*, "phantom_path: ", phantom_path   
    !infile = phantom_path // "flrw.in" 
    ! Convert file name from CCTK STRING to a fortran string for use
    ! This is required because CCTK KEYWORDS and STRINGS are passed as c pointers
    ! TODO Catch error for path longer than string length 
    call CCTK_FortranString(pathstringlength,phantom_path,path)
    infile = trim(path) // "flrw.in"
    print*, "Phantom path is: ", path 
    print*, "Infile is: ", infile
    call startrun(infile,logfile,evfile,dumpfile)
    print*, "Calling die!!"
    call die
    call CCTK_INFO("Setting pressure from initial density")
    ! Setup pressure from density using a polytrope
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)

                press(i,j,k) = poly_k * (rho(i,j,k)**poly_gamma)

            enddo 
        enddo 
    enddo 

end subroutine HelloHydro_init

subroutine calc_press(CCTK_ARGUMENTS)
    implicit none
    integer :: i,j,k
    real :: poly_k, poly_gamma
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    
    poly_k = 1.e-20
    poly_gamma = 2.0

    call CCTK_INFO("Setting pressure from density")
    ! Setup pressure from density using a polytrope
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                
                    press(i,j,k) = poly_k * (rho(i,j,k)**poly_gamma)


            enddo 
        enddo 
    enddo 



end subroutine calc_press 
