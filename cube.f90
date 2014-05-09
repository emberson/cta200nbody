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

    call dump_particles !G3

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

    subroutine initial_conditions ! Group 3
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

        ! given x,v for each particle (initial conditions)
        
        ! real, dimension(6, npmax) :: xv[ncube, ncube, *]
        
        ! from x determine the subcube of each particle
        ! also given F for each subcube     
        
        ! real, dimension(3, ngrid, ngrid, ngrid) :: force3
        
        ! update v and then x for each particle
        
!        for i=0, num_images()
     
        if (this_image() == 1) then
           write(*,*) shape(xv)
        endif

        do i = 1, npnode
           x = xv(1,i)
           y = xv(2,i)
           z = xv(3,i)

           x = mod(floor(x),ngrid)
           if (x .eq. 0) then 
              x = ncube
           end if

            y = mod(floor(y),ngrid)
           if (y .eq. 0) then 
              y = ncube
           end if

            z = mod(floor(z),ngrid)
           if (z .eq. 0) then 
              z = ncube
           end if

           xv(4,i) = xv(4,i) + force3(1,x,y,z)
           xv(1,i) = xv(1,i) + xv(4,i)

           xv(5,i) = xv(5,i) + force3(2, x,y,z)
           xv(2,i) = xv(2,i) + xv(5,i)

           xv(6,i) = xv(6,i) + force3(3, x,y,z)
           xv(3,i) = xv(3,i) + xv(6,i)
        enddo

    end subroutine update_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine send_particles !Group 2
        !
        ! Send all particles that have moved out of each node's physical volume to the appropriate neighbouring node.
        !
     
        implicit none

        integer :: i, X_new, Y_new, Z_new
        real :: x, y, z, L

        iout = 0 !number of particles that have left the subcube
        i = 1 !particle index

        L = ngrid * ncube !total size of cube

        do while (i .le. npnode) !check ith particle
           x = mod(xv(1,i), L) !find updated coordinates normalized to the cube
           y = mod(xv(2,i), L) 
           z = mod(xv(3,i), L) 
           
           xv(1,i) = x !reassign to updated coordinates
           xv(2,i) = y
           xv(3,i) = z

           X_new = 1 + (floor(x)/ngrid) !determine new subcube location
           Y_new = 1 + (floor(y)/ngrid)
           Z_new = 1 + (floor(z)/ngrid)

           if ( (X_new .eq. myCoord(1)).and.(Y_new .eq. myCoord(2)).and.(Z_new .eq. myCoord(3)) ) then
              i = i + 1 !if particle is in same subcube, check next particle
           else !if particle is not in same subcube
              write (*,*) this_image(), '#', i, 'from (',myCoord(1), ',', myCoord(2), ',', myCoord(3), ')'
              write (*,*) this_image(), '***    moved to(',X_new, ',', Y_new, ',', Z_new, ')'

              npnode[X_new, Y_new, Z_new] = npnode[X_new, Y_new, Z_new] + 1 !add the particle at its new subcube location
              xv(:,npnode[X_new, Y_new, Z_new])[X_new, Y_new, Z_new] = xv(:,i) !add particles coordinates at its new subcube location
              xv(:,i) = xv(:,npnode) !remove particle coordinates from old subcube
              npnode = npnode - 1 !remove particle from old subcube
              iout = iout + 1 !add leaving particle to count
           end if
         end do

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

