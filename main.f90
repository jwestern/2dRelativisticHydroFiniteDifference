	program main
	use M_initial
	use M_evolve
	use M_pars
	use M_derivs
        use scaling_gen
        use MKL_VSL_TYPE
        use MKL_VSL
        use MKL_DFTI
        use fourier_tools
	
	implicit none
	
!U will be conservative vars
!P will be primitive vars
	real*8, dimension(:,:,:), allocatable :: Unew, Pnew	
	real*8, dimension(:,:,:), allocatable :: Uold,Pold	
	real*8, dimension(:,:), allocatable :: x, y
 	real*8, dimension(:,:), allocatable :: vx_y, vy_x!, fx, fy
	real*8 :: dx,dy, dt, time
	integer :: ITE, i,j, Nv, Nsteps
	integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
	logical :: ltrace
	integer :: STAT,ierr
	character(len=32) :: arg
	real ( kind = 8 ), dimension(:,:), allocatable :: kx,ky,kmag,scaling,scalingsolve,scaltemp
	real*8 :: lf2
	real :: pi=3.14159265359
	logical, dimension(:,:), allocatable :: maskr,kreal,kimag

        !MKL DFT stuff - only handle, for passing around
        type(DFTI_DESCRIPTOR), POINTER :: FFTHandle
        integer :: cstrides(3), rstrides(3)

        !For VSL RNG
        real(kind=8) avg, sigma ! parameters of normal distribution
        integer :: lengths(2)

        TYPE (VSL_STREAM_STATE) :: stream
        integer(kind=4) errcode
        integer brng,method,numrands

!Read in process number from command line
	CALL getarg(1,arg)

!set for debug or not
	ltrace = .false.

!read pars
        call readpars


!set the max number of eqns for now
	Nv = 3

!allocate needed memmory
!Unew will have the updated values of the fields
!uold will have the old values of the fields
!x, y are coords 
! var 1 --> D or e (cons, prim)
! var 2 --> Sx or vx
! var 3, --> Sy or vy

!assume same extent for both coordinates
	
!allocate(Unew(Nv,N,N),Uold(Nv,N,N),x(N,N),y(N,N),Pnew(Nv,N,N),Pold(Nv,N,N), &
!     &           vx_y(N,N),vy_x(N,N))!,fy(N,N),fx(N,N))
	
	allocate(Unew(Nv,N,N),STAT=ierr)
	allocate(Uold(Nv,N,N),STAT=ierr)
	allocate(x(N,N),STAT=ierr)
	allocate(y(N,N),STAT=ierr)
	allocate(Pnew(Nv,N,N),STAT=ierr)
	allocate(Pold(Nv,N,N),STAT=ierr)
	allocate(vx_y(N,N),STAT=ierr)
	allocate(vy_x(N,N),STAT=ierr)!,fy(N,N),fx(N,N))
	allocate(kx(N-1,N-1),stat=ierr)
	allocate(ky(N-1,N-1),stat=ierr)
	allocate(kmag(N-1,N-1),stat=ierr)
	allocate(scaling(N-1,N-1),stat=ierr)
	allocate(scalingsolve(N-1,N-1),stat=ierr)
	allocate(maskr(N-1,N-1),stat=ierr)
	allocate(kreal(N-1,N-1),stat=ierr)
	allocate(kimag(N-1,N-1),stat=ierr)
	allocate(scaltemp(N-1,N-1),stat=ierr)

!Initializing rng stream
        !numrands = number of random numbers generated each time step
        !         = (1st loop size in forcingfield2)
        !                  *(2nd loop size in forcingfield2)
        numrands=(N-1)*(N-1)
        a = 0.0 !Average of gaussian distribution
        sigma  = 1.0    !Width of gaussian distribution
        brng=VSL_BRNG_MT19937  
        method=VSL_RNG_METHOD_GAUSSIAN_ICDF
        if (resume==0) then
                errcode=vslnewstream( stream, brng,  seed )
        elseif (resume==1) then  
                errcode=vslloadstreamf( stream, 'stream_state' )
        endif

        if (nskip==0) then
                write(*,*) 'RNG status: No block splitting.'
        else
                errcode=vslskipaheadstream( stream, nskip )
        endif
		
!define coords and related	
	
	dx= Lx/(N-1)
	dy = dx
	dt = cfl * dx
	
	do i = 1, N
	x(i,:) =  (i-1)*dx
	y(:,i) =  (i-1)*dy
	end do

!Define array of wavenumbers, arranged in PACK format, in grid units
!i.e. Nyquist = (N-1)/2
        call kxky(kx,ky,N-1)

!Define mask arrays which pick out real/imag parts of PACK format
        call real_imag_masks(kreal,kimag,N-1)

if (gaussforce==1) then
!define enveloping profile for Gaussian*k^2 profile
        kmag=2*pi*sqrt( kx**2 + ky**2 )/Lx
        lf2=lf*dx
        scaling=sqrt(exp(-lf2**2*kmag**2/2.))

elseif (shiftrectforce==1) then
!define enveloping profile for rectangular profile
        kmag=kx**2+ky**2
        maskr=(kmag.ge.kforcel**2).and.(kmag.lt.kforceg**2)
        scaling=1.
        where(.not.maskr) scaling=0.
endif

!Now do FFT
        lengths(1)=N-1
        lengths(2)=N-1

        STAT = DftiCreateDescriptor(FFTHandle, DFTI_DOUBLE,&
                        & DFTI_REAL, 2, lengths)
        STAT = DftiSetValue(FFTHandle, DFTI_PLACEMENT, DFTI_NOT_INPLACE)
        STAT = DftiSetValue(FFTHandle, DFTI_CONJUGATE_EVEN_STORAGE,&
                        & DFTI_COMPLEX_REAL)
        STAT = DftiSetValue(FFTHandle, DFTI_PACKED_FORMAT,&
                        & DFTI_PACK_FORMAT)
        STAT = DftiSetValue(FFTHandle, DFTI_FORWARD_SCALE, 1.0/(N-1.)**2. )
	rstrides = [0, 1, lengths(1)]
	cstrides = [0, 1, lengths(2)]

	STAT = DftiSetValue(FFTHandle, DFTI_INPUT_STRIDES, rstrides)
	STAT = DftiSetValue(FFTHandle, DFTI_OUTPUT_STRIDES, cstrides)

        STAT = DftiCommitDescriptor(FFTHandle)

!Recale 'scaling' to have desired (1,1) component in real space
        call rescale(scaling,weos**2.*F0*lf2**2.)	

!define initial data
	call initial(Uold,Pold,dx,dy,x,y,Nv,N)

if (resume==1) time = 1.0d0*nint(resume_time/dt)*dt
if (resume==0) time = 0.0d0

	Unew = Uold
	Pnew = Pold

	if(ltrace) print*, 'about to integrate'
	
	   do ITE=1, Nt
	    !integrate

        if(ltrace) print*, 'call evolve'
	    call evolve(unew,uold,Pnew,Pold,dx,dt,time,x,y,Nv,N,ITE,freq,arg,&
		&	kx,ky,scaling,errcode,method,stream,numrands,avg,sigma,&
		&	FFTHandle,lengths)	
        if(ltrace) print*, 'back from evolve'

          if (mod(ITE-1,freq).eq.0) then
          print*, 'call out', time
!calculate vorticity
	  call derivs(pold(ivy,:,:),vy_x,dx,N,1)
	  call derivs(pold(ivx,:,:),vx_y,dx,N,2)
	  	  	
          ret = gft_out_brief(trim(arg)//'vort',time, (/n,n/), 2, vy_x-vx_y)  
          ret = gft_out_brief(trim(arg)//'vx',time, (/n,n/), 2, pold(ivx,:,:))
          ret = gft_out_brief(trim(arg)//'vy',time, (/n,n/), 2, pold(ivy,:,:))
          ret = gft_out_brief(trim(arg)//'eps',time, (/n,n/), 2, pold(irho,:,:))
!          ret = gft_out_brief('Sx',time, (/n,n/), 2, uold(ivx,:,:))
!          ret = gft_out_brief('Sy',time, (/n,n/), 2, uold(ivy,:,:))
!          ret = gft_out_brief('D',time, (/n,n/), 2, uold(irho,:,:))
 
        if(ltrace) print*, 'end call out'
	  end if

	   !update the dime and save the new values
	   !in the old array. 

            time = time + dt
            Uold = Unew
	    Pold = Pnew

	   end do

          ret = gft_out_brief(trim(arg)//'vx',time, (/n,n/), 2, pold(ivx,:,:))
          ret = gft_out_brief(trim(arg)//'vy',time, (/n,n/), 2, pold(ivy,:,:))
          ret = gft_out_brief(trim(arg)//'eps',time, (/n,n/), 2, pold(irho,:,:))
!          ret = gft_out_brief('Sx',time, (/n,n/), 2, uold(ivx,:,:))
!          ret = gft_out_brief('Sy',time, (/n,n/), 2, uold(ivy,:,:))
!          ret = gft_out_brief('D',time, (/n,n/), 2, uold(irho,:,:))

	  STAT = vslsavestreamf( stream, 'stream_state' )
	  STAT = DftiFreeDescriptor(FFTHandle)

	end
