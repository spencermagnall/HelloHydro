#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine HelloHydro_init(CCTK_ARGUMENTS)
    !use io,              only:id,master,nprocs,set_io_unit_numbers,die
    !use initial,         only:initialise,finalise,startrun,endrun
    !use evolve,          only:evol_init, evol_step
    use einsteintk_wrapper
    use einsteintk_utils
    use metric_utils
    implicit none
    character(len=500) :: infile,logfile,evfile,dumpfile,path
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
    !infile = trim(path) // "flrw.in"
    infile = 'flrw.in'
    !infile = trim(infile)
    print*, "Phantom path is: ", path 
    print*, "Infile is: ", infile
    ! Use system call to copy phantom files to simulation directory
    ! This is a digusting temporary fix
    call SYSTEM('cp ~/phantom/phantomET/test/flrw* ./')
    !call set_io_unit_numbers
    !call startrun(infile,logfile,evfile,dumpfile)
    !call evol_init(infile,logfile,evfile,dumpfile)
    !call evol_step(infile,logfile,evfile,dumpfile)
    !print*, "Calling die!!"
    !call die
    ! print*, "cctk_bbox: ", cctk_bbox 
    ! print*, "cctk_nghostzones: ", cctk_nghostzones
    ! stop 
    call CCTK_INFO("Setting up metric grid in phantom")
    !print*, "Boundary size: ", boundary_size_x_lower
    call init_et2phantomgrid(cctk_gsh(1),cctk_gsh(2),cctk_gsh(3), & 
    cctk_delta_space(1), cctk_delta_space(2), cctk_delta_space(3), &
    cctk_origin_space(1),cctk_origin_space(2), cctk_origin_space(3))
    
    call CCTK_INFO("Passing the metric from ET to phantom")
    call HelloHydro_metric_global(CCTK_ARGUMENTS)
    call print_etgrid
    call CCTK_INFO("Initialising phantom and starting run")
    ! When calling phantom we also want to pass the dt
    ! of einstein toolkit as our maximun sph timestep
    call init_et2phantom(infile,CCTK_DELTA_TIME)

    !call Tmunu_init(CCTK_ARGUMENTS)
    !print*,"Tmunu ET values are: ", eTtt(6,6,6), eTxx(6,6,6), eTyy(6,6,6), eTzz(6,6,6)  
    
    call CCTK_INFO("Setting pressure from initial density")
        ! Setup pressure from density using a polytrope
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)

                rho(i,j,k) = 0.                
                !press(i,j,k) = poly_k * (rho(i,j,k)**poly_gamma)

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

subroutine Tmunu_init(CCTK_ARGUMENTS)
    !use einsteintk_wrapper
    use einsteintk_utils, only: tmunugrid
    use metric_utils
    implicit none
    integer :: i,j,k
    DECLARE_CCTK_ARGUMENTS
    DECLARE_CCTK_FUNCTIONS
    DECLARE_CCTK_PARAMETERS
    
    call CCTK_INFO("Getting initial stress energy tensor from phantom grid!")
    print*, "tmunugrid 0,0 ", tmunugrid(1,1,6:8,6:8,6:8)
    ! Add stress energy contribution to TmunuBase 
    do k=1, cctk_lsh(3)
        do j=1, cctk_lsh(2)
            do i=1, cctk_lsh(1)
                
                ! I really hate the way ET splits up its variables 
                ! This is super confusing using zero based indexing in
                ! Phantom but not in ET
                ! I Need to think of a better way to do this 
                ! As I should not directly be setting the Tmunu
                ! Scalar part of Tmunu               
                eTtt(i,j,k) = tmunugrid(1,1,i,j,k)
               
                ! Vector part of Tmunu
                eTtx(i,j,k) =  tmunugrid(1,2,i,j,k) 
                eTty(i,j,k) =  tmunugrid(1,3,i,j,k)
                eTtz(i,j,k) =  tmunugrid(1,4,i,j,k) 
                ! Tensor part ov Tmunu
                eTxx(i,j,k) =  tmunugrid(2,2,i,j,k)
                eTxy(i,j,k) =  tmunugrid(2,3,i,j,k)
                eTxz(i,j,k) =  tmunugrid(2,4,i,j,k)
                eTyy(i,j,k) =  tmunugrid(3,3,i,j,k)
                eTyz(i,j,k) =  tmunugrid(3,4,i,j,k)
                eTzz(i,j,k) =  tmunugrid(4,4,i,j,k)



            enddo 
        enddo 
    enddo 
end subroutine Tmunu_init
