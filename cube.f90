! ----------------------------------------------------------------------------------------------------
! cube.f90: Parallel N-body code written May 2014 by CTA200 students.
! Compile with ifort -fpp -c mkl_fftvec.f90 ; ifort -fpp -mkl -coarray -coarray-num-images=8 cube.f90 mkl_fftvec.o -o cube.x 
! ----------------------------------------------------------------------------------------------------

program cube

    implicit none

    !! Simulation parameters
    integer, parameter :: ngrid = 2 !each subcube has 2x2x2 mesh
    integer, parameter :: ncube = 2 !number of subcubes in one dimension
    integer, parameter :: np = 32 !number of particles
    integer, parameter :: timesteps = 3
    integer, parameter :: npmax = 4*np/ncube**3
    integer, parameter :: npen = ngrid/ncube

    !! Particle positions and velocities
    real, dimension(6, npmax) :: xv[ncube, ncube, *]!* is the number of images (8) divided by ncube twice (so will be ncube=2 again)

    !! Local density field 
    real, dimension(ngrid, ngrid, ngrid) :: rho !local to processor you are on (not a coarray)
    real, dimension(ngrid, npen, ncube, npen, ncube) :: rhol !rho redimensionalized, unpacking x and y coords
    equivalence(rho, rhol) !physical space in memory for rho and rhol the same (=> equivalent in fortran)

    !! Local force field
    real, dimension(3, ngrid, ngrid, ngrid) :: force3 !3D forcefield, components of force in x,y,z

    !! Force kernel
    complex, dimension(npen, ncube, ncube, ngrid/2+1, npen, 3) :: cforce3[ncube, ncube, *] !complex version of force3

    !! Working copies of rho for pencil routines - Group 1
    real, dimension(ngrid, npen, ncube, npen, ncube) :: rho3[ncube, ncube, *]
    real, dimension(ngrid*ncube+2, npen, ncube, npen) :: rhox[ncube, ncube, *]
    complex, dimension(ngrid*ncube/2+1, npen, ncube, npen) :: crhox[ncube, ncube, *]
    complex, dimension(npen, ncube, ncube, npen, ngrid/2+1) :: crhoy[ncube, ncube, *]
    complex, dimension(npen, ncube, ncube, ngrid/2+1, npen) :: crhoz[ncube, ncube, *]

    !! Temporary coarrays needed in the pencil routines - Group 1
    complex, dimension(ngrid/2+1, npen, ncube, npen) :: crhotmpyxglobal[ncube, ncube, *]
    complex, dimension(npen, npen, ngrid/2+1) :: crhotmpzyglobal[ncube, ncube, *]
    integer, codimension[ncube,ncube,*] :: ixr,iyr,ibr

    !! Temporary coarray needed in the Poisson solver routine - Group 1
    complex, dimension(npen, ncube, ncube, ngrid/2+1, npen) :: crhoztmp[ncube, ncube, *]

    !! Number of particles on each node and coarray to keep track of sending
    integer npnode[ncube, ncube, *] !number of particles in a given subcube
    integer iout[ncube, ncube, *] !Group 2

    !! Time step counter
    integer :: it !tracks timestep

    !! Image coordinates in the volume decomposition 
    integer mycoord(3)

    ! ----------------------------------------------------------------------------------------------------
    ! MAIN
    ! ----------------------------------------------------------------------------------------------------

    !! Only proceed if the number of images is equal to the number of cubes
    if (num_images() .ne. ncube**3) then
        write(*,*) "ERROR: num_images != ncube^3 ... ", num_images(), ncube**3
        stop
    endif

    !! Only proceed if the grid can be evenly split into pencils
    if (mod(ngrid,ncube) .ne. 0) then
        write(*,*) "ERROR: Cannot evenly decompose the mesh into pencils ... ", ngrid, ncube
        stop
    endif

    !! Store cubic image coordinate
    mycoord = this_image(xv)

    !! Start with an equal number of particles on each node
    npnode = np/ncube**3 !number of particles in subcube

    call setup_kernels
    call initial_conditions !Group 3 - randomize particles

    do it = 1, timesteps !now proceed through timesteps

        call calculate_rho !calculate density field, Group 3
        call poisson_solve !Group 1
        call update_particles !Group 2
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
        ! Store the force kernel in cforce3. 
        !
    
        implicit none

    end subroutine setup_kernels

    ! ----------------------------------------------------------------------------------------------------

    subroutine initial_conditions
        !
        ! Initialize particles to start with random positions within the range [0, ngrid*ncube] and with
        ! zero velocity, randomly fill up xv
        !

        implicit none

    end subroutine initial_conditions

    ! ----------------------------------------------------------------------------------------------------

    subroutine calculate_rho !Group 3
        !
        ! Interpolate particles onto the mesh and accumulate their density normalized to the mean density.
        !

        implicit none

    end subroutine calculate_rho

    ! ----------------------------------------------------------------------------------------------------

    subroutine poisson_solve
        !
        ! We want to solve Poisson's equation for the 3D force field. In Fourier space we have that
        ! phi(k) = -4*pi*G*rho(k)/k^2 and F(k) = -m*(k.phi(k)). So we forward transform rho(x) to 
        ! solve for rho(k) and then perform 3 inverse transforms to get the components of F(x).
        !

        implicit none

        integer :: idim

        !
        ! Fourier transform density field
        !

        rho3 = rhol
        call pencilfftforward !Group 1
        crhoztmp = crhoz

        !
        ! Get each component of F(x) in real space
        !

        do idim = 1, 3
            crhoz = crhoztmp * conjg(cforce3(:,:,:,:,:,idim))
            call pencilfftbackward !Group 1
            rhol = rho3
            force3(idim,:,:,:)=rho
        enddo

    end subroutine poisson_solve

    ! ----------------------------------------------------------------------------------------------------

    subroutine pencilfftforward !Group 1
        !
        ! Start with the real space density field rho3 and transform this to Fourier space in variable crhoz. 
        ! This involves a transformation from a cubic to pencil decomposition.
        !

        implicit none

        !
        ! First unpack the cubic representation of rho3 into a pencil representation in rhox.
        ! Pencils have their longest axis in the x dimension and their shortest in z. Pencils
        ! are stacked along y first and then z.
        !

        !! DO SOME THINGS HERE TO GO FROM RHO3 TO RHOX

        call fftvec(rhox, ngrid*ncube, ngrid**2/ncube, 1) !mkl fft

        !
        ! Now transpose pencils in the x-y plane so that they have their longest axis in y and 
        ! their shortest in z. Pencils are stacked along the z axis first and then x.
        !

        crhox = cmplx(rhox(::2,:,:,:), rhox(2::2,:,:,:))

        !! DO SOME THINGS HERE TO GO FROM CRHOX -> CRHOY

        call cfftvec(crhoy, ngrid*ncube, ngrid*(ngrid/2+1)/ncube, 1)

        !
        ! Now transpose pencils in the y-z plane so that they have their longest axis in z and shortest in y.
        ! Pencils are stacked along the x axis first and then y.
        !

        call cfftvec(crhoz, ngrid*ncube, ngrid*(ngrid/2+1)/ncube, 1)

    end subroutine pencilfftforward

    ! ----------------------------------------------------------------------------------------------------

    subroutine pencilfftbackward
        !
        ! Start with the Fourier transformed density field in crhoz and inverse transform this back to
        ! real space in variable rho3. This involves a transformation from a pencil to cubic decomposition.
        !

        implicit none

    end subroutine pencilfftbackward

    ! ----------------------------------------------------------------------------------------------------

    subroutine update_particles !Group 2
        !
        ! First update each particle's velocity by adding to it the acceleration it feels from the force 
        ! of the cell it sits in. Then move each particle according to its updated velocity and apply
        ! periodic boundary conditions.
        !

        implicit none

    end subroutine update_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine send_particles !Group 2
        !
        ! Send all particles that have moved out of each node's physical volume to the appropriate neighbouring node.
        !

      
        implicit none

    end subroutine send_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine dump_particles !Group 3
        !
        ! Write out the positions of all particles contained in the volume to a single text file. 
        ! Use only the master image for this. (processor with this_image() = 1)
        !

        implicit none

    end subroutine dump_particles

    ! ----------------------------------------------------------------------------------------------------

end program cube

