! requires dfti module, build with
! ifort -c /software/compilers/Intel/2013-sp1-14.0/composer_xe_2013_sp1.0.080/mkl/include/mkl_dfti.f90 
! then compile ifort -mkl=sequential -fpp -DSELFTEST mkl_fftvec.f90
#ifdef SELFTEST
integer, parameter :: n=4
complex z(n/2+1)
real a(n),b(n)
equivalence(a,z)
common z,b

a=(/ 1,2,3,4/)
b=a/2
write(*,*) a
call fftvec(a,n,2,1)
write(*,*) z,a
call fftvec(a,n,2,-1)
write(*,*) a
end

#endif

subroutine cfftvec(crho,nx,ny,idir)
  use mkl_dfti
  complex crho(nx*ny)
  integer Status
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
!  common /mklfft/ Desc_Handle_Dim1, Desc_Handle_r2c
  integer nx0,ny0
  logical firsttime
  save firsttime,nx0,ny0,Desc_Handle_Dim1
  data firsttime /.true./

  if (firsttime) then
     firsttime = .false.
!     write(*,*) 'init ',nx,ny
     nx0=nx
     ny0=ny
     Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_SINGLE,&
          DFTI_COMPLEX, 1, nx )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, ny )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, nx )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, nx )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1./nx )
     Status = DftiCommitDescriptor( Desc_Handle_Dim1 )  
     if (Status .ne. 0) then
        if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
           print *, 'FFT init Error: ', DftiErrorMessage(status)
           pause
        endif
     endif
  end if
  if (nx0 .ne. nx .or. ny0 .ne. ny) then
     write(*,*) 'nx,ny,nx0,ny0 inconsistency ',nx,ny,nx0,ny0 
     stop
  end if
  if (idir .eq. 1) then
     Status = DftiComputeForward( Desc_Handle_Dim1, crho )
  else
     Status = DftiComputeBackward( Desc_Handle_Dim1, crho )
  end if
  if (Status .ne. 0) then
     if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
        print *, 'dir=',idir,'c2c fft Error: ', DftiErrorMessage(status)
        pause
     endif
  endif
end subroutine cfftvec


subroutine fftvec(crho,nx,ny,idir)
  use mkl_dfti
  complex crho((nx/2+1)*ny)
  integer Status
  type(DFTI_DESCRIPTOR), POINTER :: Desc_Handle_Dim1
  integer nx0,ny0
  logical firsttime
  save firsttime,nx0,ny0,Desc_Handle_Dim1
  data firsttime /.true./

  if (firsttime) then
     firsttime = .false.
     nx0=nx
     ny0=ny
     Status = DftiCreateDescriptor(Desc_Handle_Dim1, DFTI_SINGLE,&
          DFTI_REAL, 1, nx )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_NUMBER_OF_TRANSFORMS, ny )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_INPUT_DISTANCE, nx+2 )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_OUTPUT_DISTANCE, nx+2 )
     Status = DftiSetValue( Desc_Handle_Dim1, DFTI_BACKWARD_SCALE, 1./nx )
     Status = DftiCommitDescriptor( Desc_Handle_Dim1 )     
     if (Status .ne. 0) then
        if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
           print *, 'FFT init Error: ', DftiErrorMessage(status)
           pause
        endif
     endif
  end if
  if (nx0 .ne. nx .or. ny0 .ne. ny) then
     write(*,*) 'nx,ny,nx0,ny0 inconsistency ',nx,ny,nx0,ny0 
     stop
  end if
  if (idir .eq. 1) then
     Status = DftiComputeForward( Desc_Handle_Dim1, crho )
  else
     Status = DftiComputeBackward( Desc_Handle_Dim1, crho )
  end if
  if (Status .ne. 0) then
     if (.not. DftiErrorClass(status,DFTI_NO_ERROR)) then
        print *, 'r2c fft Error: ', DftiErrorMessage(status)
        pause
     endif
  endif
end subroutine fftvec
