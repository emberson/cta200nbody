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
    integer, parameter :: timesteps = 10!3
    integer, parameter :: npmax = 4*np/ncube**3
    integer, parameter :: npen = ngrid/ncube
    integer, parameter :: nwrite = 2

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

    !integer :: temp(3),test(3) !testing index_global
    integer :: i


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
    myCoord = this_image(xv)

    !temp = (/1,2,1/)
    !test = index_global(temp)
    !write(*,*) this_image(), 'testing',test

    !! Start with an equal number of particles on each node
    npnode = np/ncube**3 !number of particles in subcube

    npnode=0
    if (this_image() == 1) npnode=1

    call setup_kernels
    call initial_conditions !Group 3 - randomize particles
    call init_xv_test

    force3(1,:,:,:)= .1
    !force3(2,:,:,:)= .05
!    force3(2,:,:,:)=.5

    do it = 1, timesteps !now proceed through timesteps

        !call calculate_rho !calculate density field, Group 3
        !call poisson_solve !Group 1
        call update_particles !Group 2
        call send_particles
	call dump_particles

        do i = 1, npnode
           write(*,*) '(x,y,z) = ', xv(1,i), xv(2,i), xv(3,i)
        end do
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

        integer :: i, j, k, idim
        real :: r2, dr3(3)

        do idim = 1, 3
            do k = 1, ngrid
                do j = 1, ngrid
                    do i = 1, ngrid
                        dr3 = delta_r(index_global((/i, j, k/)))
                        r2  = sum(dr3**2)
                        if (r2 < 0.5) then
                            rho(i, j, k) = 0.
                            cycle
                        end if
                        rho(i, j, k) = -0.1 * dr3(idim) / r2**1.5
                   enddo
                enddo
            enddo

            rho3 = rhol
            call pencilfftforward
            cforce3(:,:,:,:,:,idim) = crhoz

        enddo

    end subroutine setup_kernels

    function index_global(index_local)
        !
        ! Takes local mesh coordinates and transforms them to global mesh coordinates.
        !

        implicit none

        integer :: index_local(3), index_global(3),j,toglobal
        

        do j=1,3
           index_global(j) =(mycoord(j)-1)*ngrid+index_local(j)
        enddo

    end function index_global

    function delta_r(index_glob)
        !
        ! Takes global mesh coordinates as an input and maps this to delta_r based on the following: 
        !   * Each dimension of index_global ranges from [1, ngrid*ncube] and this is linearly mapped onto
        !     delta_r with range [0, 1, ..., ngrid*ncube/2 -1, 1.e10, -ngrid*ncube/2 + 1, ..., -2, -1]. 
        !

        implicit none

        integer :: index_glob(3)
        real delta_r(3)

        real :: L, dx, dy, dz

        L = ngrid * ncube
        dx = index_glob(1)
        dy = index_glob(2)
        dz = index_glob(3)
        
        if(dx .gt. (L - dx)) dx = dx - L
        if(dy .gt. (L - dy)) dy = dy - L
        if(dz .gt. (L - dz)) dz = dz - L
	
	if (dx .eq. L/2) dx=10**10
	if (dy .eq. L/2) dy=10**10
        if (dz .eq. L/2) dz=10**10
	
        delta_r = (/dx,dy,dz/)

    end function delta_r


    ! ----------------------------------------------------------------------------------------------------

    subroutine initial_conditions ! Group 3
        !
        ! Initialize particles to start with random positions within the range [0, ngrid*ncube] and with
        ! zero velocity, randomly fill up xv
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

    subroutine calculate_rho !Group 3
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

!        ix =  mycoord(3)
!        ixlen = ngrid/2
!        if (ix .eq. ncube) ixlen = ixlen + 1  !! Nyquist frequency (top row in z mapped to far x end here)
!        crhoy = 0.
!
!        do iy0 = 1, ncube
!            iy = mod(iy0+ix, ncube) + 1
!            ixr[iy,mycoord(1),mycoord(2)]=ix
!            sync all
!            crhotmpyxglobal= !!
!            sync all
!            crhotmpyx= !! 
!            do i=1,ixlen
!                crhoy(:,:,iy,:,i)=crhotmpyx(i,:,:,:)
!            enddo
!        enddo


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
       ! ! real space in variable rho3. This involves a transformation from a pencil to cubic decomposition.
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

        integer :: i
        real :: x,y,z
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

        ! For each particle in a node
        do i = 1, npnode
           x = xv(1,i)
           y = xv(2,i)
           z = xv(3,i)

           if ( (x.ge. (myCoord(1))*ngrid) .or. (x.lt. (myCoord(1)-1)*ngrid) .or. (y.ge. (myCoord(2))*ngrid) .or. (y.lt. (myCoord(2)-1)*ngrid) &
                      .or. (z.ge. (myCoord(3))*ngrid) .or. (z.lt. (myCoord(3)-1)*ngrid) ) then 
                      write(*,*) 'ERROR !!! Your particle is outside !!!!'
           endif

           x = mod(floor(x),ngrid)+1 !find its locations in terms of local coordinates
           y = mod(floor(y),ngrid)+1
           z = mod(floor(z),ngrid)+1

           xv(4,i) = xv(4,i) + force3(1,x,y,z) !update velocity : v(t +dt) = a(t) *dt + v(t)
           xv(1,i) = xv(1,i) + xv(4,i) !update position : x(t+dt) = x(t) + v_x(t)

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

           ! if particle is in same subcube, check next particle
           if ( (X_new .eq. myCoord(1)).and.(Y_new .eq. myCoord(2)).and.(Z_new .eq. myCoord(3)) ) then
                i = i + 1
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
         sync all !synchronize

    end subroutine send_particles

    ! ----------------------------------------------------------------------------------------------------

    subroutine dump_particles !Group 3
        !
        ! Write out the positions of all particles contained in the volume to a single text file. 
        ! Use only the master image for this. (processor with this_image() = 1)
        !

        implicit none
	integer :: i, j, k, l, digit, fileNum
	Character(LEN = 1024) :: strNum, strDigit 
		
	! Writes to a file every (nwrite)th iteration (including the first iteration).
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

    subroutine init_xv_test

      implicit none
      integer :: i

      do i = 1,npnode

         xv(1,i) = (myCoord(1)-1)*ngrid + 0.5*(1+i)
         xv(2,i) = (myCoord(2)-1)*ngrid + 0.5
         xv(3,i) = (myCoord(3)-1)*ngrid + 0.5

         end do 

      end subroutine init_xv_test
end program cube

