! compile with 
!ifort    -fpp -O3  -coarray cube.f90 -mkl mkl_fftvec.o -coarray-num-images=8
!
! ngrid is the number of grid cells on each node
integer, parameter :: ngrid=32,ncube=2,np=128*8,npmax=np*4/ncube**3
real, dimension(6,npmax) :: xv[ncube,ncube,*]
real, dimension(ngrid,ngrid,ngrid) :: rho
real, dimension(3,ngrid,ngrid,ngrid) :: force3
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



mycoord=this_image(xv)
npnode=np/ncube**3
if (num_images() .ne. ncube**3) pause 'ncube .ne. num_images'
call setup_kernels
call initial_conditions
do it=1,100
   rho=0
   call calcrho
   rho=rho-np/(1.*ngrid*ncube)**3
   call poisson_solve
   call updatexv
   call send_particles
end do
call dump_particles
contains
  function index_global(i3)
    implicit none
    integer i3(3), index_global(3)
    index_global=(mycoord-1)*ngrid+i3
  end function index_global
  function delta_r(i3)
    real delta_r(3)
    integer i3(3)
    delta_r=i3-1
    where(delta_r>ngrid*ncube/2)delta_r=delta_r-ngrid*ncube
    where(abs(delta_r)>ngrid*ncube/2.-0.5) delta_r=1.e10
  end function delta_r
  subroutine setup_kernels
    real dr3(3)
    do idim=1,3
       do k=1,ngrid
          do j=1,ngrid
             do i=1,ngrid
                dr3=delta_r(index_global((/i,j,k/)))
                r2=sum(dr3**2)
                if (r2<0.5) then
                   rho(i,j,k)=0
                   cycle
                end if
                rho(i,j,k)=-0.1*dr3(idim)/r2**1.5
             end do
          end do
       end do
       rho3=rhol
!       call gather_matrix
       do inode=1,8*0
          sync all
          if (idim .eq. 1 .and. this_image() .eq. inode) then
             write(*,*) inode,mycoord
             write(*,'(4(4f6.3,1x))') rho
          end if
       end do
       call pencilfftforward(rho3,rhox,crhox,crhoy,crhoz,ngrid,ncube)
!       write(*,*) idim,this_image(),maxval(abs(crhoz)),maxval(rho3)
       cforce3(:,:,:,:,:,idim)=crhoz
    end do
  end subroutine setup_kernels
subroutine initial_conditions
  integer, parameter :: nseedmax=200
  integer iseed(nseedmax)

  call random_seed()
  call random_seed(size=iseedsize)
  if (size>nseedmax) pause 'need bigger nseedmax'
  call random_seed(get=iseed(:iseedsize))
  iseed=iseed+this_image()
  call random_seed(put=iseed(:iseedsize))
  xv(:,:npnode)=0
  call random_number(xv(:3,:npnode))
  xv(:3,:npnode)=(xv(:3,:npnode)+spread(mycoord-1,2,npnode))*ngrid+1
end subroutine initial_conditions
subroutine calcrho
  do i=1,npnode
     ix3=xv(:3,i)-(mycoord-1)*ngrid
     if (any( floor((ix3-1.) / ngrid) .ne. 0)) then
        write(*,*) 'out of bound',ix3,mycoord
        cycle
     end if
     rho(ix3(1),ix3(2),ix3(3))= rho(ix3(1),ix3(2),ix3(3))+1   
  end do
  sum1=sum(rho)
end subroutine calcrho
subroutine updatexv
  integer ix3(3)
!  call random_number(xv(4:,:npnode))
!  xv(4:,:npnode)=2*xv(4:,:npnode)-1

  do ipart=1,npnode
     ix3=xv(:3,ipart)-(mycoord-1)*ngrid
     xv(4:,ipart)=xv(4:,ipart)+force3(:,ix3(1),ix3(2),ix3(3))
  end do
  xv(:3,:npnode)=xv(:3,:npnode)+xv(4:,:npnode)
  xv(:3,:npnode)=modulo(xv(:3,:npnode)+ngrid*ncube-1,ngrid*ncube*1.)+1
  if (this_image() .eq. 1 .and. mod(it,10) .eq. 1) write(*,*) sum1-np/ncube**3,maxval(xv(4,:npnode)),minval(xv(4,:npnode))!,maxval(force3(1,:,:,:)),minval(force3(1,:,:,:))

end subroutine updatexv
integer function islab(x)
  real x
  islab=floor((x-1)/ngrid)+1
end function islab
subroutine send_particles
  integer itarget3(3)
  do idim=1,3
     do idir=-1,1,2
        iout=npmax+1
        ipart=1
        do 
           if (ipart > npnode) exit
           if ( any( (islab(xv(idim,ipart))-mycoord(idim))*idir &
                .eq. (/1,1-ncube/)) ) then
!              write(*,*) idim,idir,xv(:,ipart),mycoord
              iout=iout-1
              xv(:,iout)=xv(:,ipart)
              xv(:,ipart)=xv(:,npnode)
              npnode=npnode-1
           else
              ipart=ipart+1
           end if
        end do
        sync all
        itarget3=mycoord
        itarget=mod(mycoord(idim)+idir+ncube-1,ncube)+1
        itarget3(idim)=itarget
        nptarget=npnode[itarget3(1),itarget3(2),itarget3(3)]+1
        npnewtarget=nptarget+npmax-iout
        nplim=2*npmax-iout[itarget3(1),itarget3(2),itarget3(3)]
        if (npnewtarget .gt. nplim) pause 'too many particles'
        xv(:,nptarget:npnewtarget)[itarget3(1),itarget3(2),itarget3(3)]=&
             xv(:,iout:npmax)
        npnode[itarget3(1),itarget3(2),itarget3(3)]=npnewtarget
        sync all    
     end do
  end do
end subroutine send_particles

subroutine poisson_solve
!  rho=0
!  if (this_image() .eq. 1) rho(1,1,1)=1
  rho3=rhol
  call pencilfftforward(rho3,rhox,crhox,crhoy,crhoz,ngrid,ncube)
  crhoztmp=crhoz
  do inode=1,8*0
     sync all
     if (this_image() .eq. inode) write(*,*) this_image(),crhoy
  end do
  do idim=1,3
     crhoz=crhoztmp*conjg(cforce3(:,:,:,:,:,idim))
     call pencilfftbackward(rho3,rhox,crhox,crhoy,crhoz,ngrid,ncube)
     rhol=rho3
!     call gather_matrix
     force3(idim,:,:,:)=rho
  end do
!  write(*,*) this_image(),maxval(force3)
end subroutine poisson_solve

subroutine pencilfftforward(rho3,rhox,crhox,crhoy,crhoz,ngrid,ncube)
! think of every matrix as a collection of ngrid**3*ncube**3 object
! index always positive cyclic (ijk) (abc) (xyz)
! input matrix rho3(i,a,j,b,k,c)[x,y,z]
! rhox(i,a,x,j,b,k)[y,c,z] but needs 2 extra, so
! rhox(i*a*x+2,j,b,k)[]
! for y,z keep i*a together
! complex crhoy(j,b,y,k,i*a/2+1)[c,z,x]
! complex crhoz(k,c,z,i*a/2+1,j)[x,b,y]
!
! we transpose the relevant orderings
  real, dimension(ngrid/ncube*ncube,ngrid/ncube,ncube,ngrid/ncube,ncube) :: rho3[ncube,ncube,*]
  real, dimension(ngrid/ncube*ncube*ncube+2,ngrid/ncube,ncube,ngrid/ncube) :: rhox[ncube,ncube,*]
  complex, dimension(ngrid*ncube/2+1,ngrid/ncube,ncube,ngrid/ncube) :: crhox[ncube,ncube,*]
  complex, dimension(ngrid/ncube,ncube,ncube,ngrid/ncube,ngrid/2+1) :: crhoy[ncube,ncube,*]
  complex, dimension(ngrid/ncube,ncube,ncube,ngrid/2+1,ngrid/ncube) :: crhoz[ncube,ncube,*]
  complex, dimension(ngrid/2+1,ngrid/ncube,ncube,ngrid/ncube) :: crhotmpyx
  complex, dimension(ngrid/ncube,ngrid/ncube,ngrid/2+1) :: crhotmpzy
  integer mycoord(3)

  mycoord=this_image(rho3)
  if (mod(ngrid,ncube).ne.0) pause 'ERROR: ncube does not divide ngrid'
! align x-axis first
  ic=mycoord(2)
  iy=mycoord(1)
  iz=mycoord(3)
!  write(*,*) 'start first transpose',rhox
  rhox=0
  sync all
  do ix0=1,ncube ! for unknown reasons, ifort is very fast on this loop
! even without synchronization within each iteration
     ix=mod(ix0+ic,ncube)+1  ! stagger the communications
     rhox((ix-1)*ngrid+1:ix*ngrid,:,:,:)=rho3(:,:,:,:,ic)[ix,iy,iz]
  end do
!  write(*,*) 'done first transpose'
  call fftvec(rhox,ngrid*ncube,ngrid**2/ncube,1)
  crhox=cmplx(rhox(::2,:,:,:),rhox(2::2,:,:,:))
! now transpose x-y
  ix=mycoord(3)
  ixlen=ngrid/2
  if (ix .eq. ncube) ixlen=ixlen+1  ! nyquist frequency
  crhoy=0
!  sync all
  do iy0=1,ncube
     iy=mod(iy0+ix,ncube)+1
     ixr[iy,mycoord(1),mycoord(2)]=ix
     sync all
     crhotmpyxglobal=crhox((ixr-1)*ngrid/2+1:ixr*ngrid/2+1,:,:,:)
     sync all
     crhotmpyx=crhotmpyxglobal[iy,mycoord(1),mycoord(2)]
     do i=1,ixlen
        crhoy(:,:,iy,:,i)=crhotmpyx(i,:,:,:)
     end do
  end do
  call cfftvec(crhoy,ngrid*ncube,ngrid*(ngrid/2+1)/ncube,1)
! now transpose z
  ix=mycoord(1)
  ib=mycoord(2)
  iy=mycoord(3)
  sync all
  do ic0=1,ncube
     ic=mod(ic0+ib,ncube)+1
     do iz0=1,ncube
        iz=mod(iz0+iy,ncube)+1
        ibr[ic,iz,ix]=ib
        iyr[ic,iz,ix]=iy
        sync all
        crhotmpzyglobal=crhoy(:,ibr,iyr,:,:)![ic,iz,ix]
        sync all
        crhotmpzy=crhotmpzyglobal[ic,iz,ix]        
        do j=1,ngrid/ncube
           crhoz(:,ic,iz,:,j)=crhotmpzy(j,:,:)
        end do
     end do
  end do
  call cfftvec(crhoz,ngrid*ncube,ngrid*(ngrid/2+1)/ncube,1)
end subroutine pencilfftforward
subroutine pencilfftbackward(rho3,rhox,crhox,crhoy,crhoz,ngrid,ncube)
  real, dimension(ngrid/ncube*ncube,ngrid/ncube,ncube,ngrid/ncube,ncube) :: rho3[ncube,ncube,*]
  real, dimension(ngrid/ncube*ncube*ncube+2,ngrid/ncube,ncube,ngrid/ncube) :: rhox[ncube,ncube,*]
  complex, dimension(ngrid*ncube/2+1,ngrid/ncube,ncube,ngrid/ncube) :: crhox[ncube,ncube,*]
  complex, dimension(ngrid/ncube,ncube,ncube,ngrid/ncube,ngrid/2+1) :: crhoy[ncube,ncube,*]
  complex, dimension(ngrid/ncube,ncube,ncube,ngrid/2+1,ngrid/ncube) :: crhoz[ncube,ncube,*]
  complex, dimension(ngrid/2+1,ngrid/ncube,ncube,ngrid/ncube) :: crhotmpyx
  complex, dimension(ngrid/ncube,ngrid/ncube,ngrid/2+1) :: crhotmpzy
  integer mycoord(3)

  mycoord=this_image(rho3)

  ix=mycoord(1)
  ib=mycoord(2)
  iy=mycoord(3)
  sync all
  call cfftvec(crhoz,ngrid*ncube,ngrid*(ngrid/2+1)/ncube,-1)

  do ic0=1,ncube
     ic=mod(ic0+ib,ncube)+1
     do iz0=1,ncube
        iz=mod(iz0+iy,ncube)+1
        do j=1,ngrid/ncube
           crhotmpzy(j,:,:)=crhoz(:,ic,iz,:,j)
        end do
        crhotmpzyglobal[ic,iz,ix]=crhotmpzy
        ibr[ic,iz,ix]=ib
        iyr[ic,iz,ix]=iy
        sync all
        crhoy(:,ibr,iyr,:,:)=crhotmpzyglobal
     end do
  end do

  call cfftvec(crhoy,ngrid*ncube,ngrid*(ngrid/2+1)/ncube,-1)
! now transpose x-y
  ix=mycoord(3)

  do iy0=1,ncube
     iy=mod(iy0+ix,ncube)+1
     ixlen=ngrid/2
     if (ix .eq. ncube) ixlen=ixlen+1  ! nyquist frequency
     do i=1,ixlen
        crhotmpyx(i,:,:,:)=crhoy(:,:,iy,:,i)
     end do
     crhotmpyxglobal[iy,mycoord(1),mycoord(2)]=crhotmpyx
     ixr[iy,mycoord(1),mycoord(2)]=ix
     sync all
     ixlen=ngrid/2
     if (ixr .eq. ncube) ixlen=ixlen+1  ! nyquist frequency
     imin=(ixr-1)*ngrid/2+1
     imax=imin+ixlen-1
     crhox(imin:imax,:,:,:)=crhotmpyxglobal(:ixlen,:,:,:)
  end do
  call fftvec(crhox,ngrid*ncube,ngrid**2/ncube,-1)
  rhox(::2,:,:,:)=real(crhox)
  rhox(2::2,:,:,:)=aimag(crhox)

  ic=mycoord(2)
  iy=mycoord(1)
  iz=mycoord(3)
  do ix=1,ncube
     rho3(:,:,:,:,ic)[ix,iy,iz]=rhox((ix-1)*ngrid+1:ix*ngrid,:,:,:)
  end do
  sync all
end subroutine pencilfftbackward

subroutine dump_particles()
  if (this_image() .ne. 1) return
     open(10,file='x.dat')
     do i=1,ncube
        do j=1,ncube
           do k=1,ncube
              write(10,'(3f10.2)') xv(:3,:npnode[i,j,k])[i,j,k]
           end do
        end do
     end do     
end subroutine dump_particles

  subroutine gather_matrix
    sync all
    if (this_image() .ne. 1) return
    write(*,*) ' '
    do k=1,ncube
       do j=1,ncube
          do i=1,ncube
             rhol=rho3[i,j,k]
             rhoglobal(:,i,:,j,:,k)=rho
          end do
       end do
    end do
!    write(*,*) sum(rhoglobalflat)
    write(*,'(4(4f6.3,1x))') rhoglobalflat
!    open(10,file='transpose.dat')
!    write(10,'(4(4f6.0,5x))') rhoglobal
  end subroutine gather_matrix

end

#ifdef NR2
subroutine fftvec(rho1,nx,ny,idir)
  call fftvec2(rho1,rho1,nx,ny,idir)
end subroutine fftvec
subroutine fftvec2(rho1,crho,nx,ny,idir)
  complex crho(nx/2+1,ny)
  real rho1(nx+2,ny)
  complex crho2(nx)
  do j=1,ny
     if (idir .eq. 1) then
        crho2=0
        crho2=rho1(:nx,j)
        call four1(crho2,nx,idir)
        crho(:,j)=crho2(:nx/2+1)
     else
        crho2(:nx/2+1)=crho(:,j)
        crho2(nx/2+2:)=conjg(crho2(nx/2:2:-1))
        call four1(crho2,nx,idir)
        rho1(:nx,j)=real(crho2)/nx
     end if
  end do
end subroutine fftvec2

subroutine cfftvec(crho,nx,ny,idir)
  complex crho(nx,ny)
  do j=1,ny
     call four1(crho(1,j),nx,idir)
  end do
  if (idir .eq. -1) crho=crho/nx
end subroutine cfftvec
#endif
