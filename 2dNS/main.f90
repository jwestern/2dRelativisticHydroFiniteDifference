	program main
	use M_initial
	use M_evolve
	use M_pars
	use M_derivs
	use POISSON_SOLVER
	use scaling_gen
	use MKL_VSL_TYPE
	use MKL_VSL
	use MKL_DFTI
	use fourier_tools

	implicit none
	
!U will be conservative vars
!P will be primitive vars
	real*8, dimension(:,:,:), allocatable :: Unew, Uold
!	real*8, dimension(:,:,:), allocatable :: Uold,Pold	
	real*8, dimension(:,:), allocatable :: x, y
 	real*8, dimension(:,:), allocatable :: vx, vy, psi, vort
	real*8 :: dx,dy, dt, time
	integer :: ITE, i,j, Nv, Nsteps
        integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
	logical :: ltrace
	integer :: STAT,ierr
	character(len=32) :: arg
	real ( kind = 8 ), dimension(:,:), allocatable :: kx,ky,kmag,scaling,msk,scalingsolve,scaltemp
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

	real*8, dimension(:,:), allocatable :: tn, dxT, dxdxT, dyT, dydyT

	real*8 :: alph
!Read in process number from command line
	CALL getarg(1,arg)

!set for debug or not
	ltrace = .false.

!read pars
        call readpars


!set the max number of eqns for now
	Nv = 1

!allocate needed memmory
!Unew will have the updated values of the fields
!uold will have the old values of the fields
!x, y are coords 
! var 1 --> D or e (cons, prim)
! var 2 --> Sx or vx
! var 3, --> Sy or vy

 	allocate(Unew(Nv,N,N),STAT=ierr)
    allocate(Uold(Nv,N,N),STAT=ierr)
    allocate(x(N,N),STAT=ierr)
    allocate(y(N,N),STAT=ierr)
!    allocate(Pnew(Nv,N,N),STAT=ierr)
!    allocate(Pold(Nv,N,N),STAT=ierr)
    allocate(vx(N,N),STAT=ierr)
    allocate(vy(N,N),STAT=ierr)!,fy(N,N),fx(N,N))
    allocate(psi(N,N),STAT=ierr)
    allocate(vort(N,N),STAT=ierr)
!    allocate(vy2(N,N),STAT=ierr)
!    allocate(vx2(N,N),STAT=ierr)
!    allocate(dxvy2(N,N),STAT=ierr)
!    allocate(dyvx2(N,N),STAT=ierr)
    allocate(kx(N-1,N-1),stat=ierr)
    allocate(ky(N-1,N-1),stat=ierr)
    allocate(kmag(N-1,N-1),stat=ierr)
    allocate(scaling(N-1,N-1),stat=ierr)
	allocate(scalingsolve(N-1,N-1),stat=ierr)
    allocate(maskr(N-1,N-1),stat=ierr)
	allocate(msk(N-1,N-1),stat=ierr)
	allocate(kreal(N-1,N-1),stat=ierr)
	allocate(kimag(N-1,N-1),stat=ierr)
	allocate(scaltemp(N-1,N-1),stat=ierr)
	allocate(tn(N,N),stat=ierr)
	allocate(dxT(N,N),stat=ierr)
	allocate(dxdxT(N,N),stat=ierr)
	allocate(dyT(N,N),stat=ierr)
	allocate(dydyT(N,N),stat=ierr)

!Initializing rng stream
	!numrands = number of random numbers generated each time step
	!         = (1st loop size in forcingfield2)
	!		   *(2nd loop size in forcingfield2)
	numrands=(N-1)*(N-1)
	a = 0.0	!Average of gaussian distribution
	sigma  = 1.0	!Width of gaussian distribution
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

!Define the two-thirds filtering mask
	msk=1.
	do i=1,N-1
		do j=1,N-1
			if (sqrt(kx(i,j)**2.+ky(i,j)**2.)>=float((N-1)/3)) then
!			if (abs(kx(i,j))>float((N-1)/3)) then
				msk(i,j)=0.
			endif
!			if (abs(ky(i,j))>float((N-1)/3)) then
!				msk(i,j)=0.
!			endif
		enddo
	enddo
!write(*,*) kx
!write(*,*) ky

if (gaussforce==1) then
!define enveloping profile for Gaussian*k^2 profile
        kmag=2*pi*sqrt( kx**2 + ky**2 )/Lx
        lf2=lf*dx
        scaling=sqrt(exp(-lf2**2*kmag**2/2.))*kmag**2.

elseif (bofforce==1) then
!define enveloping profile for Gaussian profile
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

!Apply 2/3-dealiasing filter to 'scaling'. Note it's already in Fourier space
	scaling=msk*scaling

!Recale 'scaling' to have desired (1,1) component in real space
	call rescale(scaling,F0)
!!!!!!!!!!!!TESTING
!	tester=0.
!	tn=sin(2*pi*1*(x+y)/10.) !+ &
!	&sin(2*pi*15*(x)/10.)*sin(2*pi*15*(y)/10.)
!	alph=0.5
!	tn=(1/alph**4.)*((x-5.)**2.+(y-5.)**2.-2.*alph**2.)&
!	&*exp(-0.5*((x-5.)/alph)**2-0.5*((y-5.)/alph)**2)
!	call derivs(tn,dxT,dx,N,1)
!	call derivs(dxT,dxdxT,dx,N,1)
!	call derivs(tn,dyT,dx,N,2)
!	call derivs(dyT,dydyT,dx,N,2)
	
!	call derivs(dxdxT,dxT,dx,N,1)
!	call derivs(dydyT,dyT,dx,N,2)
!	call derivs(dxT,dxdxT,dx,N,1)
!	call derivs(dyT,dydyT,dx,N,2)
!	ret = gft_out_brief('d4tn',1.0d0,(/n,n/),2,dxdxT+dydyT)

!	ret = gft_out_brief('laptn',1.0d0,(/n,n/),2,dxdxT+dydyT)
!	tn=dxdxT+dydyT
!	call derivs(tn,dxT,dx,N,1)
!	call derivs(dxT,dxdxT,dx,N,1)
!	call derivs(tn,dyT,dx,N,2)
!	call derivs(dyT,dydyT,dx,N,2)
!	ret = gft_out_brief('lap2tn',1.0d0,(/n,n/),2,dxdxT+dydyT)
!	tn=sin(2*pi*1*(x+y)/10.)
!	call dissip(tn,dxT,dx,N)
!	ret = gft_out_brief('disp',1.0d0,(/n,n/),2,dxT)
!	ret = gft_out_brief('fdxt',1.0d0, (/n,n/), 2, dxT)
!	ret = gft_out_brief('t',1.0d0, (/n,n/), 2, tn)

!Determine and print the energy injection wavenumber
!XXX: THIS WORKS ONYL IF FORCING VORTICITY DIRECTLY, SO CHANGE IT/REMOVE IT
	scaltemp=scaling**2.
	where (kimag) scaltemp=0.
	Stat = DftiComputeBackward(FFTHandle, scaltemp(:,1), scalingsolve(:,1) )
	scalingsolve=scalingsolve*F0/scalingsolve(1,1)
	tn(1:n-1,1:n-1)=scalingsolve
	tn(n,:)=tn(1,:)
	tn(:,n)=tn(:,1)
	call inverse_laplacian(tn,N,dx,kx,ky,FFTHandle,lengths)
	write(*,*) 1./sqrt(maxval(tn)/F0)
!
!
!
!
!define initial data
	call initial(Uold,dx,dy,x,y,Nv,N)

if (resume==1) time = 1.0d0*nint(resume_time/dt)*dt+dt
if (resume==0) time = 0.0d0
	Unew = Uold

	if(ltrace) print*, 'about to integrate'
	
	   do ITE=1, Nt
	    !integrate

        if(ltrace) print*, 'call evolve'

	    call evolve(unew,uold,dx,dt,time,x,y,Nv,N,ITE,freq,arg,kx,ky,scaling,&
					&errcode,method,stream,numrands,avg,sigma,&
					&FFTHandle,lengths,msk)

        if(ltrace) print*, 'back from evolve'

          if (mod(ITE-1,freq).eq.0) then
          print*, 'call out', time

          ret = gft_out_brief(trim(arg)//'vort',time, (/n,n/), 2, unew(ivort,:,:))
	  vort=unew(ivort,:,:)

	  call inverse_laplacian(vort,N,dx,kx,ky,FFTHandle,lengths) !now it becomes stream function
	  call derivs(vort,vx,dx,N,2) !get partial_y psi
	  call derivs(-vort,vy,dx,N,1) !get -partial_x psi

          ret = gft_out_brief(trim(arg)//'vx',time, (/n,n/), 2, vx)
          ret = gft_out_brief(trim(arg)//'vy',time, (/n,n/), 2, vy)

        if(ltrace) print*, 'end call out'
	  end if

	   !update the time and save the new values
	   !in the old array. 

            time = time + dt
            Uold = Unew

	   end do
		STAT = vslsavestreamf( stream, 'stream_state' )
		STAT = DftiFreeDescriptor(FFTHandle)

!         ret = gft_out_brief(trim(arg)//'vort',time, (/n,n/), 2, unew(ivort,:,:))

!	  call inverse_laplacian(vort,N,dx,kx,ky) !now it becomes stream function
!	  call derivs(vort,vx,dx,N,2) !get partial_y psi
!	  call derivs(-vort,vy,dx,N,1) !get -partial_x psi
	  
!          ret = gft_out_brief(trim(arg)//'vx',time, (/n,n/), 2, vx)
!          ret = gft_out_brief(trim(arg)//'vy',time, (/n,n/), 2, vy)

	end
