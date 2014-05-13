! ----------------------------------------------------------------------------------------------------
! cube.f90: Parallel N-body code written May 2014 by CTA200 students.
! Compile with ifort -fpp -c mkl_fftvec.f90 ; ifort -fpp -mkl -coarray -coarray-num-images=8 cube.f90 mkl_fftvec.o -o cube.x 
! ----------------------------------------------------------------------------------------------------

program cube

    implicit none

    !! Simulation parameters
    integer, parameter :: ngrid = 4
    integer, parameter :: ncube = 2
    integer, parameter :: np = 32
    integer, parameter :: timesteps = 3
    integer, parameter :: npmax = 4*np/ncube**3
    integer, parameter :: npen = ngrid/ncube

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

    complex, dimension(ngrid*ncube/2+1, npen, ncube, npen) :: crhox_copy[ncube, ncube, *]


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

    do it = 1, 1!timesteps

!        call calculate_rho
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

    end subroutine initial_conditions

    ! ----------------------------------------------------------------------------------------------------

    subroutine calculate_rho
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
        
        !
        ! Initialize rho3
        !


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
        integer i,j,k, jStart,jStop,IMG,ind
        complex, dimension(ngrid/2+1, npen, ncube, npen) :: temp
        !
        ! First unpack the cubic representation of rho3 into a pencil representation in rhox.
        ! Pencils have their longest axis in the x dimension and their shortest in z. Pencils
        ! are stacked along y first and then z.
        !
        IMG=this_image()
        
        ind = (IMG-1)*ngrid**3 + 1
        do i=1,ngrid
           do  j=1,ngrid
              do k=1,ngrid
                 rho(k,j,i)=ind
                 ind = ind + 1!nin-1)*ncube**3+(i-1)*ngrid**2+(j-1)*ngrid+k
              enddo
           enddo
        enddo
        rho3=rhol
        
        !sync all
        !call sleep(IMG)
        !write(*,*) 'IMG rho3 = ',IMG,' rho3 =' ,rho3
        !call sleep(num_images()-IMG)

        ! GO FROM RHO3 TO RHOX
        rhox = -IMG
        do i = 1,ncube
           rhox((i-1)*ngrid+1:i*ngrid,:,:,:)=rho3(:,:,:,:,mycoord(2))[i,mycoord(1),mycoord(3)]
        enddo
        sync all
        call sleep(IMG)
        write(*,*)  'IMG rhox = ',IMG,' rhox = ', rhox
        call sleep(num_images()-IMG)
       
       ! call fftvec(rhox, ngrid*ncube, ngrid**2/ncube, 1)

        !
        ! Now transpose pencils in the x-y plane so that they have their longest axis in y and 
        ! their shortest in z. Pencils are stacked along the z axis first and then x.
        !
        
        crhox = cmplx(rhox(::2,:,:,:), rhox(2::2,:,:,:))
 
!        call sleep(IMG)       
!        write(*,*)  'IMG crhox* = ',IMG,' crhox* = ', crhox
!        call sleep(num_images()-IMG)
!        sync all
        

        !complex, dimension(ngrid*ncube/2+1, npen, ncube, npen) :: crhox[ncube, ncube, *]
        !complex, dimension(npen, ncube, ncube, npen, ngrid/2+1) :: crhoy[ncube, ncube, *]
        !complex, dimension(npen, ncube, ncube, ngrid/2+1, npen) :: crhoz[ncube, ncube, *]

        ! GO FROM CRHOX -> CRHOY
        crhoy=0
        !jStart=(mycoord(3)-1)*ngrid/2+1
        jStart=1
        !jStop=ngrid*mycoord(3)/2
        jStop=ngrid/2
        !write(*,*)jStart, ' <--jStart  jStop--> ',jStop
        if (mycoord(3) == ncube) jStop= jStop+1
        do i = 1,ncube
           temp=crhox((mycoord(3)-1)*ngrid/2+1:mycoord(3)*ngrid/2+1,:,:,:)[i,mycoord(1),mycoord(2)]
           ! temp=crhox(1:ngrid,:,:,:)[i,mycoord(1),mycoord(2)]
           do j = jStart, jStop
              crhoy(:,:,i,:,j)=temp(j,:,:,:)!crhox(j,:,:,:)[i,mycoord(1),mycoord(2)]
           enddo
        enddo
!        sync all
!        call sleep(IMG)
!        write(*,*) 'IMG crhoy = ',IMG,' crhoy = ', crhoy
!        call sleep(num_images()-IMG)
        !call cfftvec(crhoy, ngrid*ncube, ngrid*(ngrid/2+1)/ncube, 1)

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
        integer i,j,k, jStart,jStop,IMG
        complex, dimension(npen,ncube,npen,ngrid/2+1) :: temp

        call cfftvec(crhoz, ngrid*ncube, ngrid*(ngrid/2+1)/ncube, -1)

!        crhox_copy(:,:,:,:) = crhox(:,:,:,:)
        !
        ! Transpose from crhoz to crhoy so that pencils have their longest axis in y and 
        ! their shortest in z. Pencils are stacked along the z axis first and then x.
        !

        !! PUT CODE HERE FOR CRHOZ -> CRHOY

        !call cfftvec(crhoy, ngrid*ncube, ngrid*(ngrid/2+1)/ncube, -1)

        !
        ! Transpose from crhoy to crhox so that pencils have their longest axis in x and 
        ! their shortest in z. Pencils are stacked along y first and then z.
        !
    
        
        !! PUT CODE HERE FOR CRHOY -> CRHOX
        IMG=this_image()
        jStart=1
        jStop=ngrid/2
        !write(*,*)jStart, ' <--jStart  jStop--> ',jStop
        !if (mycoord(3) == ncube) jStop= jStop+1
        do i = 1,ncube
           !temp=crhoy(:,:,:,:,(mycoord(3)-1)*ngrid/2+1:mycoord(3)*ngrid/2+1)[mycoord(2),mycoord(3),i]
           
           !do j = jStart, jStop
            !  crhox(j,:,:,:)=temp(:,:,i,:,j)!crhox(j,:,:,:)[i,mycoord(1),mycoord(2)]
              
           !enddo
           temp = crhoy(:,:,myCoord(1),:,:)[mycoord(2),mycoord(3),i]
           if (i == ncube) jStop = jStop +1
           do j = jStart,jStop

              crhox((i-1)*ngrid/2+j,:,:,:) = temp(:,:,:,j)
           enddo
        enddo

        sync all
!        crhox = crhox - crhox_copy

        sync all
        call sleep(IMG)
        write(*,*)  'IMG crhox = ',IMG,' crhox = ', crhox
        !call sleep(num_images()-IMG)


        !call fftvec(crhox, ngrid*ncube, ngrid**2/ncube, -1)

        rhox(::2,:,:,:) = real(crhox)
        rhox(2::2,:,:,:) = aimag(crhox)

        !
        ! Here unpack the pencil decomposition in rhox to the cubic decomposition in rho3.
        !

        !! PUT CODE HERE FOR RHOX -> RHO3

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

    end subroutine dump_particles

    ! ----------------------------------------------------------------------------------------------------

end program cube

