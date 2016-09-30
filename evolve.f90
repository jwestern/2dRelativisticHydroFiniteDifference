	module M_EVOLVE
	use MKL_VSL_TYPE
	use MKL_VSL
	implicit none
	
	contains
	
	
	subroutine evolve(unew,uold,Pnew,Pold,dx,dt,time,x,y,Nv,N,ITE,freq,arg,kx,ky,scaling,&
			&errcode,method,stream,numrands,avg,sigma,FFTHandle,lengths)
	use m_derivs
	use m_boundary
	use m_rhs
	use m_pars, only : STMETH, disip, irho, ivx,ivy,weos,printforcing
	use fourier_tools
	use MKL_DFTI
	implicit none
	
	real*8, dimension(:,:,:), intent(inout):: unew,Pnew
	real*8, dimension(:,:,:), intent(in):: uold,Pold
	real*8, dimension(:,:), intent(in) :: kx,ky,scaling
	integer, intent(in) :: ITE,freq
	type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
	integer, dimension(:), intent(inout) :: lengths
	character(len=*), intent(in) :: arg

!VSL RNG stuff
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
	real*8, dimension(:,:,:), allocatable :: Fxu,FyU,dxFxU,dyFyU
	real*8, dimension(:,:), allocatable :: forcingx,forcingy,streamfunc
	real*8 EPS, factor2,factor3
	integer i,j, gft_out_brief, ret
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
	&	 Fxu(Nv,N,N),FyU(Nv,N,N),dxFxU(Nv,N,N),dyFyU(Nv,N,N),&
	&	 forcingx(N,N),forcingy(N,N),streamfunc(N,N))

!Create stream function for forcing field	
	call forcingfield2(streamfunc,N,scaling,errcode,method,stream,numrands,avg,sigma,FFTHandle,lengths)

!Define the forcing field in the usual way, as the "curl" of the stream function (\partial_y,-\partialx)\psi
	call derivs( streamfunc,forcingx,dx,N,2)
	call derivs(-streamfunc,forcingy,dx,N,1)

        if (mod(ITE-1,freq).eq.0) then
        if(printforcing.eq.1)then
	write(*,*) 'printing force'
        ret = gft_out_brief(trim(arg)//'fx',time, (/n,n/), 2, forcingx)
        ret = gft_out_brief(trim(arg)//'fy',time, (/n,n/), 2, forcingy)
        endif
        endif

! We'll do the rk 3rd stuff. The procedure is basically
! the same 3 times, (i) evaluate derivs; (ii) evaluate rhs;
! (iii) get the new value
	
	EPS = disip 

!!!!!!!!!!!!!!!!!!! FIRST STEP	
        call floor(uold(irho,:,:),N)
	
!get the fluxes
	call fluxes(uold,Pold,FxU,FyU,N,Nv)
!calculate derivs	
	do i = 1,Nv
	  call derivs(uold(i,:,:),dxU(i,:,:),dx,N,1)
	  call derivs(uold(i,:,:),dyU(i,:,:),dx,N,2)
	  call derivs(FxU(i,:,:),dxFxU(i,:,:),dx,N,1)
	  call derivs(FyU(i,:,:),dyFyU(i,:,:),dx,N,2)
	  call dissip(uold(i,:,:),disu(i,:,:),dx,N)
	end do

	printforcing=0
	call rhs(uold,dxU,dyU,dxFxU,dyFyU,urk1,dx,N,&
	&	pold,FFTHandle,lengths)
	
        if(ltrace) print*, 'call boun'
	call boundary(uold,dxU, dyU, urk1,dx,time,x)
	
	urk1 = urk1 + EPS*disu
	
	unew = uold + dt * urk1
	unew(ivx,:,:) = unew(ivx,:,:) + sqrtdt*forcingx
	unew(ivy,:,:) = unew(ivy,:,:) + sqrtdt*forcingy
	
	call floor(unew(irho,:,:),N)
	do i=1,N
	 do j=1,N
           call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
	 end do
	end do

!!!!!!!!!!!!!!!!!!! SECOND STEP	

!get the fluxes
	call fluxes(unew,Pnew,FxU,FyU,N,Nv)
!calculate derivs	
	do i = 1,Nv
	  call derivs(unew(i,:,:),dxU(i,:,:),dx,N,1)
	  call derivs(unew(i,:,:),dyU(i,:,:),dx,N,2)
	  call derivs(FxU(i,:,:),dxFxU(i,:,:),dx,N,1)
	  call derivs(FyU(i,:,:),dyFyU(i,:,:),dx,N,2)
	  call dissip(unew(i,:,:),disu(i,:,:),dx,N)	  
	end do
!	print *, 'made it up to 2nd rhs call'
	printforcing=0
	call rhs(unew,dxU,dyU,dxFxU,dyFyU,urk2,dx,N,&
	&	pold,FFTHandle,lengths)
!	print *, 'made it past 2nd rhs call'
	call boundary(unew,dxU, dyU, urk2,dx,time,x)
	
	urk2 = urk2 + EPS*disu

	unew = uold + 0.5 * dt * ( urk1 + urk2 )
	unew(ivx,:,:) = unew(ivx,:,:) + sqrtdt*forcingx
	unew(ivy,:,:) = unew(ivy,:,:) + sqrtdt*forcingy
	
	call floor(unew(irho,:,:),N)
	do i=1,N
	 do j=1,N
           call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
	 end do
	end do
!!!!!!!!!!!!!!!!!!! THIRD STEP	

!get the fluxes
!	call fluxes(unew,Pnew,FxU,FyU,N,Nv)
!calculate derivs	
!	do i = 1,Nv
!	  call derivs(unew(i,:,:),dxU(i,:,:),dx,N,1)
!	  call derivs(unew(i,:,:),dyU(i,:,:),dx,N,2)
!	  call derivs(FxU(i,:,:),dxFxU(i,:,:),dx,N,1)
!	  call derivs(FyU(i,:,:),dyFyU(i,:,:),dx,N,2)
!	  call dissip(unew(i,:,:),disu(i,:,:),dx,N)	  
!	end do
	
!	printforcing=1
!	call rhs(unew,dxU, dyU,dxFxU,dyFyU, urk3,dx,N,x,pold,time,dt,printforcing,ITE,freq,arg)
!print *, 'made it past 3rd rhs call'
	
!	call boundary(unew,dxU, dyU, urk3,dx,time,x)
	
!	urk3 = urk3 + EPS*disu

!	unew = uold + dt * urk3
	
!	call floor(unew(irho,:,:),N)
!	do i=1,N
!	 do j=1,N
!          call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
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
!	  call derivs(FxU(i,:,:),dxFxU(i,:,:),dx,N,1)
!	  call derivs(FyU(i,:,:),dyFyU(i,:,:),dx,N,2)
!	  call dissip(unew(i,:,:),disu(i,:,:),dx,N)	  
!	end do

!	call rhs(unew,dxU, dyU,dxFxU,dyFyU, urk4,dx,N,x,pold,time,dt,printforcing,ITE,freq,arg)

!	call boundary(unew,dxU, dyU, urk4,dx,time,x)

!	urk4 = urk4 + EPS*disu

!        unew = uold + dt/6. * (urk1+2.*urk2+2.*urk3+urk4)
!   else 	
!	unew = uold + dt/9. * (2.*urk1+3.*urk2+4.*urk3)		
!   end if

!	call floor(unew(irho,:,:),N)
!	do i=1,N
!	 do j=1,N
!          call contoprim(Unew(:,i,j),Pnew(:,i,j),i,j,weos)
!	 end do
!	end do

!deallocate memmory
	deallocate(dxU,dyU,urk1,urk2,urk3,urk4,disu, &
	&	    Fxu,FyU,dxFxU,dyFyU,forcingx, &
	&	    forcingy,streamfunc)

	
	end subroutine evolve
	
	
!FLUXES HERE FOR SIMPLICITY
	subroutine fluxes(u,P, Fx,Fy,N,Nv)
	use m_pars, only : irho,ivx,ivy
	
	implicit none
	real*8, dimension(:,:,:) :: Fx, Fy, U, P
	integer N, Nv
	real*8, dimension(:,:), pointer :: eps,vx,vy,D,Sx,Sy
	target :: U,P
	
	
!define the mapping
	D => u(irho,:,:)
	Sx => u(ivx,:,:)
	Sy => u(ivy,:,:)	
		
	eps => P(irho,:,:)
	vx => P(ivx,:,:)
	vy => P(ivy,:,:)

!flux for rho eqn	
	Fx(irho,:,:) = Sx
	Fy(irho,:,:) = Sy

!flux for S_i eqn
	Fx(ivx,:,:) = eps+Sx*vx
	Fy(ivx,:,:) = Sx*vy

	Fx(ivy,:,:) = Sy*vx
	Fy(ivy,:,:) = eps+Sy*vy
		
	
	end subroutine fluxes	
	
	
	subroutine floor(u,N)
        use m_pars, only : valuefloor
	implicit none
	real*8, dimension(:,:) :: u
	integer N, i, j
	
	do i=1,N
	 do j=1,N
	 if(u(i,j).lt.valuefloor)  u(i,j) = valuefloor
	 end do
	end do

	
	end subroutine floor

	subroutine contoprim(U,P,i,j,weos)
        use m_pars, only : valuefloor
	implicit none
	real*8 :: U(3),P(3),weos
	integer :: i,j

	real*8 :: WW2, ss, sqrinWW
	
	ss = (U(2)**2 + U(3)**2)/U(1)**2
	
	sqrinWW = (weos+1)**2.-4.*weos*ss
	if(ss.ge.0) then
	  sqrinWW = sqrt(sqrinWW)
	else
	  print*,'trouble, need something',i,j,ss,U(1)
		stop
	end if

	if(ss.lt.valuefloor**2) then
	   WW2 = 1.
	else
	   WW2 = (1+weos)*2*(ss-1.)/(2.*ss -weos -1. - sqrinWW)
	end if

	if(WW2.lt.0) then
		print*,'oops WW2 negative', 'WW2=', WW2, 'ss=', ss, 'U(1)=', U(1)**2, 'U(2)=', U(2)**2, 'U(3)=', U(3)**2
		stop
	endif

!now compute
	P(1) = U(1)/(-1.+(weos+1.)/WW2)
	P(2) = U(2)*WW2/((weos+1.)*P(1))
	P(3) = U(3)*WW2/((weos+1.)*P(1))

	end subroutine contoprim



	end module M_EVOLVE
