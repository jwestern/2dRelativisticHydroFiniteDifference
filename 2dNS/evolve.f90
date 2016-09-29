	module M_EVOLVE
	use MKL_VSL_TYPE
	use MKL_VSL
	implicit none
	
	contains
	
	
	subroutine evolve(unew,uold,dx,dt,time,x,y,Nv,N,ITE,freq,arg,kx,ky,scaling,&
			&errcode,method,stream,numrands,avg,sigma,FFTHandle,lengths,&
			&msk)
	use m_derivs
	use m_boundary
	use m_rhs
	use m_pars, only : STMETH,disip,ivort,weos,bofforce,gaussforce,&
			  &shiftgaussforce,shiftcosforce,F0,lf,printforcing
	use POISSON_SOLVER
	use fourier_tools
	use MKL_DFTI
	implicit none
	
	real*8, dimension(:,:,:), intent(inout):: unew
	real*8, dimension(:,:,:), intent(inout):: uold
	real*8, dimension(:,:), intent(in) :: kx,ky,scaling
	real*8 , dimension(:,:), intent(inout) :: msk
	integer, intent(in) :: ITE,freq
	type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
	integer, dimension(:), intent(inout) :: lengths

!VSL RNG stuff
	character(len=*), intent(in) :: arg
    real(kind=8), intent(in) :: avg, sigma ! parameters of normal distribution
    TYPE (VSL_STREAM_STATE), intent(inout) :: stream
    integer(kind=4), intent(inout) :: errcode
    integer, intent(in) :: method,numrands
!
	real*8, dimension(:,:):: x,y	
	real*8 dx, dt, time, sqrtdt
	integer :: Nv,N
	
!keep in some variables the rhs
	real*8, dimension(:,:,:), allocatable :: dxU,dyU,urk1,urk2,urk3,urk4,disu
	real*8, dimension(:,:,:), allocatable :: psi, dxPsi, dyPsi, lapvort, lapvortx, lapvorty
	real*8, dimension(:,:), allocatable :: forcing,streamfunc
	real*8 EPS, factor2,factor3
	integer :: i,j, gft_out_brief, ret
	logical :: ltrace

	sqrtdt=sqrt(dt)
!for debug
	ltrace = .false.
	
!allocate memmory needed to take all steps
!dxU and dyU will keep the derivatives of the field
!urk1,2,3,4 are to keep the intermidiate rhs values for the RK update
!FxU and FyU are the fluxes in the x and y direction of each U variable
!dxFxU and dyFyU are the derivatives of the x/y flux wrt the x/y coord
! we add a routine fluxes at the end to compute the fluxes
! the derivatives, rhs, boundaries have their on files where the routines
! are kept.

	allocate(dxU(Nv,N,N),dyU(Nv,N,N),urk1(Nv,N,N),urk2(Nv,N,N),&
	&        urk3(Nv,N,N),urk4(Nv,N,N),disu(Nv,N,N), &
	&	 psi(Nv,N,N),dxPsi(Nv,N,N),dyPsi(Nv,N,N),lapvort(Nv,N,N),&
	&	 lapvortx(Nv,N,N),lapvorty(Nv,N,N),forcing(N,N),streamfunc(N,N))

	call forcingfield2(streamfunc,N,scaling,errcode,method,stream,numrands,avg,sigma,FFTHandle,lengths)
	forcing=streamfunc

if (printforcing==1) then
	ret = gft_out_brief('f',time, (/n,n/), 2, forcing)
endif

! We'll do the rk 3rd stuff. The procedure is basically
! the same 3 times, (i) evaluate derivs; (ii) evaluate rhs;
! (iii) get the new value
	
	EPS = disip 

	if(STMETH.eq.3) then
	factor2 = 0.75
	else if(STMETH.eq.4) then
	factor2 = 0.5
	else
	print*, 'nocoded!'
	end if	
!!!!!!!!!!!!!!!!!!! FIRST STEP	
!        call floor(uold(ivort,:,:),N)
	
!get the fluxes
!	call fluxes(uold,Pold,FxU,FyU,N,Nv)
!calculate derivs	

	do i = 1,Nv
	  call derivs(uold(i,:,:),dxU(i,:,:),dx,N,1)
	  call derivs(uold(i,:,:),dyU(i,:,:),dx,N,2)
	  psi(i,:,:)=uold(i,:,:) !Will be replaced with psi in the next line
!write(*,*) psi(i,1,1:9)
	  call inverse_laplacian(psi(i,:,:),N,dx,kx,ky,FFTHandle,lengths)
!write(*,*) psi(i,1,1:9)
	  call derivs(psi(i,:,:),dxpsi(i,:,:),dx,N,1)
	  call derivs(psi(i,:,:),dypsi(i,:,:),dx,N,2)
	  call derivs(dxU(i,:,:),lapvortx(i,:,:),dx,N,1)
	  call derivs(dyU(i,:,:),lapvorty(i,:,:),dx,N,2)
	  call fftfilter(dxU(i,:,:),N,msk,FFTHandle,lengths)
	  call fftfilter(dyU(i,:,:),N,msk,FFTHandle,lengths)
	  call fftfilter(dxpsi(i,:,:),N,msk,FFTHandle,lengths)
	  call fftfilter(dypsi(i,:,:),N,msk,FFTHandle,lengths)
	  lapvort(i,:,:)=lapvortx(i,:,:)+lapvorty(i,:,:)
	  call dissip(uold(i,:,:),disu(i,:,:),dx,N)
	end do

	call rhs(uold,dxU, dyU,dxpsi,dypsi,lapvort,urk1,dx,N,x,&
	&	time,dt,printforcing,ITE,freq,arg,FFTHandle,lengths,msk)
	
        if(ltrace) print*, 'call boun'
!	call boundary(uold,dxU, dyU, urk1,dx,time,x)
	
	urk1 = urk1 + EPS*disu

	unew(ivort,:,:) = uold(ivort,:,:) + dt*urk1(ivort,:,:) + sqrtdt*forcing
	
!	call floor(unew(irho,:,:),N)
!	do i=1,N
!	 do j=1,N
!          call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
!	 end do
!	end do

!!!!!!!!!!!!!!!!!!! SECOND STEP	

!get the fluxes
!	call fluxes(unew,Pnew,FxU,FyU,N,Nv)
!calculate derivs	
	do i = 1,Nv
	  call derivs(unew(i,:,:),dxU(i,:,:),dx,N,1)
	  call derivs(unew(i,:,:),dyU(i,:,:),dx,N,2)
	  psi(i,:,:)=unew(i,:,:) !Will be replaced with psi in the next line
!write(*,*) psi(i,1,1:9)
	  call inverse_laplacian(psi(i,:,:),N,dx,kx,ky,FFTHandle,lengths)
!write(*,*) psi(i,1,1:9)
	  call derivs(psi(i,:,:),dxpsi(i,:,:),dx,N,1)
	  call derivs(psi(i,:,:),dypsi(i,:,:),dx,N,2)
	  call derivs(dxU(i,:,:),lapvortx(i,:,:),dx,N,1)
	  call derivs(dyU(i,:,:),lapvorty(i,:,:),dx,N,2)
          call fftfilter(dxU(i,:,:),N,msk,FFTHandle,lengths)
          call fftfilter(dyU(i,:,:),N,msk,FFTHandle,lengths)
          call fftfilter(dxpsi(i,:,:),N,msk,FFTHandle,lengths)
          call fftfilter(dypsi(i,:,:),N,msk,FFTHandle,lengths)
          lapvort(i,:,:)=lapvortx(i,:,:)+lapvorty(i,:,:)
	  call dissip(unew(i,:,:),disu(i,:,:),dx,N)	  
	end do
!	print *, 'made it up to 2nd rhs call'

	call rhs(unew,dxU, dyU,dxpsi,dypsi,lapvort,urk2,dx,N,x,&
	&	time,dt,printforcing,ITE,freq,arg,FFTHandle,lengths,msk)
!	print *, 'made it past 2nd rhs call'
!	call boundary(unew,dxU, dyU, urk2,dx,time,x)
!write(*,*) urk2(1,1,1:9)	
	urk2 = urk2 + EPS*disu
!write(*,*) urk2(1,1,1:9)
	unew(ivort,:,:) = uold(ivort,:,:) + 0.5*dt*(urk1(ivort,:,:)+urk2(ivort,:,:)) + sqrtdt*forcing
	call fftfilter(unew(ivort,:,:),N,msk,FFTHandle,lengths)
!write(*,*) unew(1,1,1:9)	
!	call floor(unew(irho,:,:),N)
!	do i=1,N
!	 do j=1,N
!           call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
!	 end do
!	end do
!!!!!!!!!!!!!!!!!!! THIRD STEP	

!get the fluxes
!	call fluxes(unew,Pnew,FxU,FyU,N,Nv)
!calculate derivs	
!	do i = 1,Nv
!	  call derivs(unew(i,:,:),dxU(i,:,:),dx,N,1)
!	  call derivs(unew(i,:,:),dyU(i,:,:),dx,N,2)
!	  psi(i,:,:)=-unew(i,:,:) !Will be replaced with psi in the next line
!write(*,*) psi(i,1,1:9)
!	  call inverse_laplacian(psi(i,:,:),N,dx)
!write(*,*) psi(i,1,1:9)
!	  call derivs(psi(i,:,:),dxpsi(i,:,:),dx,N,1)
!	  call derivs(psi(i,:,:),dypsi(i,:,:),dx,N,2)
!	  call derivs(dxU(i,:,:),lapvortx(i,:,:),dx,N,1)
!	  call derivs(dyU(i,:,:),lapvorty(i,:,:),dx,N,2)
!	  lapvort(i,:,:)=lapvortx(i,:,:)+lapvorty(i,:,:)
!	  call dissip(unew(i,:,:),disu(i,:,:),dx,N)  
!	end do
	
!	printforcing=1
!	call rhs(unew,dxU, dyU,dxpsi,dypsi,lapvort,urk3,dx,N,x,&
!	&			time,dt,printforcing,ITE,freq,arg)
!print *, 'made it past 3rd rhs call'
	
!	call boundary(unew,dxU, dyU, urk3,dx,time,x)
	
!	urk3 = urk3 + EPS*disu

!	unew = uold + dt * urk3
	
!	call floor(unew(irho,:,:),N)
!	do i=1,N
!	 do j=1,N
!           call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
!	 end do
!	end do	
!!!!!!!!!!!!!!!!!!! FOURTH STEP	

!   if (STMETH.eq.4) THEN	
!get the fluxes
!	call fluxes(unew,Pnew,FxU,FyU,N,Nv)
!calculate derivs	
!	do i = 1,Nv
!	  call derivs(unew(i,:,:),dxU(i,:,:),dx,N,1)
!	  call derivs(unew(i,:,:),dyU(i,:,:),dx,N,2)
!	  psi(i,:,:)=-unew(i,:,:) !Will be replaced with psi in the next line
!	  call inverse_laplacian(psi(i,:,:),N,dx)
!	  call derivs(psi(i,:,:),dxpsi(i,:,:),dx,N,1)
!	  call derivs(psi(i,:,:),dypsi(i,:,:),dx,N,2)
!	  call derivs(dxU(i,:,:),lapvortx(i,:,:),dx,N,1)
!	  call derivs(dyU(i,:,:),lapvorty(i,:,:),dx,N,2)
!	  lapvort(i,:,:)=lapvortx(i,:,:)+lapvorty(i,:,:)
!	  call dissip(unew(i,:,:),disu(i,:,:),dx,N)	  
!	end do

!	call rhs(unew,dxU, dyU,dxpsi,dypsi,lapvort,urk1,dx,N,x,&
!	&			time,dt,printforcing,ITE,freq,arg)

!	call boundary(unew,dxU, dyU, urk4,dx,time,x)

!	urk4 = urk4 + EPS*disu

!        unew = uold + dt/6. * (urk1+2.*urk2+2.*urk3+urk4)
!   else 	
!	unew = uold + dt/9. * (2.*urk1+3.*urk2+4.*urk3)		
!   end if

!	call floor(unew(irho,:,:),N)
!	do i=1,N
!	 do j=1,N
!           call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
!	 end do
!	end do

!deallocate memmory
	deallocate(dxU,dyU,urk1,urk2,&
	&        urk3,urk4,disu, &
	&	 psi,dxPsi,dyPsi,lapvort,&
	&	 lapvortx,lapvorty)

	
	end subroutine evolve

	subroutine init_random_seed()
	INTEGER :: i, n, clock
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed

	CALL RANDOM_SEED(size=n)
	ALLOCATE(seed(n))

	CALL SYSTEM_CLOCK(COUNT=clock)

	seed=clock+37*(/ (i-1, i=1,n) /)
	CALL RANDOM_SEED(PUT=seed)

	DEALLOCATE(seed)
	end subroutine
	
!FLUXES HERE FOR SIMPLICITY
!	subroutine fluxes(u,P, Fx,Fy,N,Nv)
!	use m_pars, only : irho,ivx,ivy
	
!	implicit none
!	real*8, dimension(:,:,:) :: Fx, Fy, U, P
!	integer N, Nv
!	real*8, dimension(:,:), pointer :: eps,vx,vy,D,Sx,Sy
!	target :: U,P
	
	
!define the mapping
!	D => u(irho,:,:)
!	Sx => u(ivx,:,:)
!	Sy => u(ivy,:,:)	
		
!	eps => P(irho,:,:)
!	vx => P(ivx,:,:)
!	vy => P(ivy,:,:)

!flux for rho eqn	
!	Fx(irho,:,:) = Sx
!	Fy(irho,:,:) = Sy

!flux for S_i eqn
!	Fx(ivx,:,:) = eps+Sx*vx
!	Fy(ivx,:,:) = Sx*vy

!	Fx(ivy,:,:) = Sy*vx
!	Fy(ivy,:,:) = eps+Sy*vy
		
	
!	end subroutine fluxes	
	
	
!	subroutine floor(u,N)
!        use m_pars, only : valuefloor
!	implicit none
!	real*8, dimension(:,:) :: u
!	integer N, i, j
	
!	do i=1,N
!	 do j=1,N
!	 if(u(i,j).lt.valuefloor)  u(i,j) = valuefloor
!	 end do
!	end do

	
!	end subroutine floor

!	subroutine contoprim(U,P,i,j,weos)
!        use m_pars, only : valuefloor
!	implicit none
!	real*8 :: U(3),P(3),weos
!	integer :: i,j

!	real*8 :: WW2, ss, sqrinWW
	
!	ss = (U(2)**2 + U(3)**2)/U(1)**2
	
!	sqrinWW = (weos+1)**2.-4.*weos*ss
!	if(ss.ge.0) then
!	  sqrinWW = sqrt(sqrinWW)
!	else
!	  print*,'trouble, need something',i,j,ss,U(1)
!		stop
!	end if

!	if(ss.lt.valuefloor**2) then
!	   WW2 = 1.
!	else
!	   WW2 = (1+weos)*2*(ss-1.)/(2.*ss -weos -1. - sqrinWW)
!	end if

!	if(WW2.lt.0) then
!		print*,'oops WW2 negative', 'WW2=', WW2, 'ss=', ss, 'U(1)=', U(1)**2, 'U(2)=', U(2)**2, 'U(3)=', U(3)**2
!		stop
!	endif

!now compute
!	P(1) = U(1)/(-1.+(weos+1.)/WW2)
!	P(2) = U(2)*WW2/((weos+1.)*P(1))
!	P(3) = U(3)*WW2/((weos+1.)*P(1))

!	end subroutine contoprim



	end module M_EVOLVE
