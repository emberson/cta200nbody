! ----------------------------------------------------------------------------------------------------
! cube.f90: Parallel N-body code written May 2014 by CTA200 students.
! Compile with 
! ----------------------------------------------------------------------------------------------------

program cube

    implicit none
   
    !! Simulation parameters
    integer, parameter :: ngrid = 2
    integer, parameter :: ncube = 2
    integer, parameter :: np = 32
    integer, parameter :: timesteps = 2
    integer, parameter :: npmax = 4*np/ncube**3

    !! Particle positions and velocities
    real, dimension(6, npmax) :: xv[ncube, ncube, *]

    !! Local density field 
    real, dimension(ngrid, ngrid, ngrid) :: rho

    !! Local force field
    real, dimension(3, ngrid, ngrid, ngrid) :: force3

    !! Misc
    integer :: it
    real sum1[ncube,ncube,*]
    integer npnode[ncube,ncube,*],iout[ncube,ncube,*]
    integer iseed(4)
    integer mycoord(3),ix3(3)
    real, dimension(ngrid,ngrid/ncube,ncube,ngrid/ncube,ncube) :: rhol
    equivalence(rho,rhol)
    real, dimension(ngrid,ncube,ngrid,ncube,ngrid,ncube) :: rhoglobal
    real, dimension(ngrid*ncube,ngrid*ncube,ngrid*ncube) :: rhoglobalflat
    equivalence(rhoglobal,rhoglobalflat)
    character*1 fn1
    real, dimension(ngrid/ncube*ncube,ngrid/ncube,ncube,ngrid/ncube,ncube) :: rho3[ncube,ncube,*]
    real, dimension(ngrid/ncube*ncube*ncube+2,ngrid/ncube,ncube,ngrid/ncube) :: rhox[ncube,ncube,*]
    complex, dimension(ngrid*ncube/2+1,ngrid/ncube,ncube,ngrid/ncube) :: crhox[ncube,ncube,*]
    complex, dimension(ngrid/ncube,ncube,ncube,ngrid/ncube,ngrid/2+1) :: crhoy[ncube,ncube,*]
    complex, dimension(ngrid/ncube,ncube,ncube,ngrid/2+1,ngrid/ncube), codimension[ncube,ncube,*] :: crhoz,crhoz2,crhoztmp
    complex, dimension(ngrid/ncube,ncube,ncube,ngrid/2+1,ngrid/ncube,3),codimension[ncube,ncube,*] :: cforce3
    complex, dimension(ngrid/2+1,ngrid/ncube,ncube,ngrid/ncube),codimension[ncube,ncube,*] :: crhotmpyxglobal
    complex, dimension(ngrid/ncube,ngrid/ncube,ngrid/2+1),codimension[ncube,ncube,*] :: crhotmpzyglobal
    integer, codimension[ncube,ncube,*] :: ixr,iyr,ibr

    ! ----------------------------------------------------------------------------------------------------
    ! MAIN
    ! ----------------------------------------------------------------------------------------------------

    !! Only proceed if the number of images is equal to the number of cubes
    if (num_images() .ne. ncube**3) then
        write(*,*) "ERROR: num_images != ncube^3 ... ", num_images(), ncube**3
        stop
    endif

    !! Store cubic image coordinate
    mycoord = this_image(xv)

    call setup_kernels
    call initial_conditions

    do it = 1, timesteps

        call calculate_rho
        call poisson_solve
        call update_particles
        call send_particles

    enddo

    call dump_particles

    ! ----------------------------------------------------------------------------------------------------
    ! SUBROUTINES
    ! ----------------------------------------------------------------------------------------------------

contains

    ! ----------------------------------------------------------------------------------------------------

    subroutine setup_kernels
        !
        ! 
        !
    
        implicit none

    end subroutine setup_kernels

    ! ----------------------------------------------------------------------------------------------------

    subroutine initial_conditions
        !
        !
        !

        implicit none

    end subroutine initial_conditions

    ! ----------------------------------------------------------------------------------------------------

    subroutine calculate_rho
        !
        !
        !

        implicit none

    end subroutine calculate_rho

    ! ----------------------------------------------------------------------------------------------------

    subroutine poisson_solve
        !
        !
        !

        implicit none

    end subroutine poisson_solve

    ! ----------------------------------------------------------------------------------------------------

    subroutine update_particles
        !
        !
        !

        implicit none

    end subroutine update_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine send_particles
        !
        !
        !

        implicit none

    end subroutine send_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine dump_particles
        !
        !
        !

        implicit none

    end subroutine dump_particles

    ! ----------------------------------------------------------------------------------------------------

end program cube

