! ----------------------------------------------------------------------------------------------------
! cube.f90: Parallel N-body code written May 2014 by CTA200 students.
! Compile with ifort -fpp -c mkl_fftvec.f90 ; ifort -fpp -mkl -coarray -coarray-num-images=8 cube.f90 mkl_fftvec.o -o cube.x 
! ----------------------------------------------------------------------------------------------------

program cube

    implicit none

    !! Simulation parameters
    integer, parameter :: ngrid = 2
    integer, parameter :: ncube = 2
    integer, parameter :: np = 32
    integer, parameter :: timesteps = 3
    integer, parameter :: npmax = 4*np/ncube**3
    integer, parameter :: npen = ngrid/ncube
    integer, parameter :: nwrite = 2

    !! Particle positions and velocities
    real, dimension(6, npmax) :: xv[ncube, ncube, *]

    !! Local density field 
    real, dimension(ngrid, ngrid, ngrid) :: rho
    real, dimension(ngrid, npen, ncube, npen, ncube) :: rhol
    equivalence(rho, rhol)

    !! Local force field
    real, dimension(3, ngrid, ngrid, ngrid) :: force3

    !! Force kernel
    complex, dimension(npen, ncube, ncube, ngrid/2+1, npen, 3) :: cforce3[ncube, ncube, *]

    !! Working copies of rho for pencil routines 
    real, dimension(ngrid, npen, ncube, npen, ncube) :: rho3[ncube, ncube, *]
    real, dimension(ngrid*ncube+2, npen, ncube, npen) :: rhox[ncube, ncube, *]
    complex, dimension(ngrid*ncube/2+1, npen, ncube, npen) :: crhox[ncube, ncube, *]
    complex, dimension(npen, ncube, ncube, npen, ngrid/2+1) :: crhoy[ncube, ncube, *]
    complex, dimension(npen, ncube, ncube, ngrid/2+1, npen) :: crhoz[ncube, ncube, *]

    !! Temporary coarrays needed in the pencil routines
    complex, dimension(ngrid/2+1, npen, ncube, npen) :: crhotmpyxglobal[ncube, ncube, *]
    complex, dimension(npen, npen, ngrid/2+1) :: crhotmpzyglobal[ncube, ncube, *]
    integer, codimension[ncube,ncube,*] :: ixr,iyr,ibr

    !! Temporary coarray needed in the Poisson solver routine
    complex, dimension(npen, ncube, ncube, ngrid/2+1, npen) :: crhoztmp[ncube, ncube, *]

    !! Number of particles on each node and coarray to keep track of sending
    integer npnode[ncube, ncube, *]
    integer iout[ncube, ncube, *]

    !! Time step counter
    integer :: it

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
    npnode = np/ncube**3

    call setup_kernels
    call initial_conditions

    do it = 1, timesteps

        call calculate_rho
        call poisson_solve
        call update_particles
        call send_particles
	call dump_particles

    enddo

    

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
        ! zero velocity.
        !

        implicit none
	integer :: i, j, seedSize, time
	integer, allocatable :: seed(:)

	! Sets all values to 0 (which will automatically handle setting all velocity components to 0).
	xv = 0
	
	! Prepares a seed to be used for generating random values.
	seedSize = -3
	call random_seed(size = seedSize)
        allocate(seed(seedSize))
	
	! Records the time, upon which the seed values for random number generation will be based.
	call system_clock(time)
	
	! Allocates seed values for random number generation, multipled by this_image() to give a different seed for each processor.
        seed(:) = time * this_image()
        call random_seed(put = seed)
		
	! All nodes iterate over the particles they contain and initialize random positions.
	do i = 1, npnode
		! Initializes random xyz position components, but only with values between 0 and 1.
		call RANDOM_NUMBER (xv(1:3,i))
	
		do j = 1, 3
			! Shifts the xyz components that were generated to the proper sub-cube this processor manages.
			xv(j,i) = (ngrid * xv(j,i)) + (ngrid * (mycoord(j) - 1))
		enddo	
	enddo

    end subroutine initial_conditions

    ! ----------------------------------------------------------------------------------------------------

    subroutine calculate_rho
        !
        ! Interpolate particles onto the mesh and accumulate their density normalized to the mean density.
        !

        implicit none

	! Make empty array for the mesh coordinates.
	integer :: mesh_coord(3), i, j

	! Resets the density to 0 every time the method is called.
	rho = 0

	! Changes the coordinates of each particle (based on the overall cube) to be in terms of the mesh of the sub-cube.
	do i = 1, npnode
		do j = 1, 3
			! This formula is essentially the reversal of the equation in initial_conditions.
			mesh_coord(j) = xv(j,i) - (ngrid * (mycoord(j) - 1))
		enddo
	
	! Increase the density of the particles within that mesh by one (particle). 
	rho(mesh_coord(1) + 1, mesh_coord(2) + 1, mesh_coord(3) + 1) = rho(mesh_coord(1) + 1, mesh_coord(2) + 1, mesh_coord(3) + 1) + 1
	enddo

	! Reduce the density of each mesh to a fraction of the total number of particles.
	rho = rho/np
	
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
        call pencilfftforward
        crhoztmp = crhoz

        !
        ! Get each component of F(x) in real space
        !

        do idim = 1, 3
            crhoz = crhoztmp * conjg(cforce3(:,:,:,:,:,idim))
            call pencilfftbackward
            rhol = rho3
            force3(idim,:,:,:)=rho
        enddo

    end subroutine poisson_solve

    ! ----------------------------------------------------------------------------------------------------

    subroutine pencilfftforward
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

        call fftvec(rhox, ngrid*ncube, ngrid**2/ncube, 1)

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

    subroutine update_particles
        !
        ! First update each particle's velocity by adding to it the acceleration it feels from the force 
        ! of the cell it sits in. Then move each particle according to its updated velocity and apply
        ! periodic boundary conditions.
        !

        implicit none

    end subroutine update_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine send_particles
        !
        ! Send all particles that have moved out of each node's physical volume to the appropriate neighbouring node.
        !

        implicit none

    end subroutine send_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine dump_particles
        !
        ! Write out the positions of all particles contained in the volume to a single text file. 
        ! Use only the master image for this.
        !

        implicit none
	integer :: i, j, k, l, digit, fileNum
	Character(LEN = 1024) :: strNum, strDigit 
		
	! Writes to a file every (nwrite)th iteration.
	if (mod(it - 1, nwrite) == 0) then 
		
		! Only the master node will write the file.
		if (this_image() == 1) then 
			
			! Writes the file number based on nwrite.
			fileNum = (it - 1) / nwrite + 1

			! Determines the number of digits in the file number, used for string formatting.
			digit = log10(real(fileNum)) + 1

			! Converts the file number to a string for file naming, and opens the file.
			write(strDigit, "(i1)") digit
			write(strNum, "(i" // trim(strDigit) // ")") fileNum
			open (unit = 1, file = "Particle" // trim(strNum) // ".txt")                             

			! Iterates over all processors and particles, and records the information (with proper column formatting).
			do i=1, ncube
				do j=1, ncube
					do k=1, ncube
						do l=1, npnode[i,j,k]
							write(1,'(6f10.2)') xv(:6,l)[i,j,k]
						enddo 
					enddo
				enddo
	
			enddo
	
			! Saves the file.
			close(1)	
		endif
	endif
	
    end subroutine dump_particles


    ! ----------------------------------------------------------------------------------------------------

end program cube

