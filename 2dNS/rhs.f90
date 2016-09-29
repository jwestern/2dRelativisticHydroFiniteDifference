 	MODULE M_RHS
	implicit none
	
	contains
	
	subroutine rhs(u,dxu,dyU,dxpsi,dypsi,lapvort,source,dx,Nx,x,time,dt,printforcing,ITE,freq,arg,&
			FFTHandle,lengths,msk)
	use m_pars, only : c, s,k_r, k_e, gamma,ivort,solenoidalforce,rhoforce,fric,relativisticfriction, &
			& forcestr, fricstr, shiftrectforce, kforcel, kforceg, kfric, filterfric,lf,F0, &
			& gaussforce,shift,width, &
			& shiftgaussforce, shiftcosforce, shift2, width2, bofforce, nu
	use fourier_tools
!	use fourier_repackage
	use M_derivs

	implicit none
	real*8, dimension(:,:,:), intent(in) :: u, dxu,dyu,dxpsi,dypsi,lapvort
	real*8, dimension(:,:,:), intent(inout) :: source
	real*8, dimension(:,:), intent(in) :: x
	real*8, intent(in) :: time, dt
	integer, intent(in) :: printforcing,ITE,freq
	character(len=*), intent(in) :: arg
	real*8, dimension(:,:), intent(inout) :: msk
	type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
        integer, dimension(:), intent(inout) :: lengths

	real*8 :: dx, v, ppi, e, wavenumber, alpha, beta, propfact
	integer :: i, j, Nx, ret, gft_out_brief
	logical :: ltrace

!My changes
	real*8, dimension(:,:), allocatable :: coordsx,coordsy,forcing,streamfunc,friction,nonlinterm
	real*8 :: extremevals(4)
	real(kind=8), dimension(:,:), allocatable :: storage,storage_x,storage_y

	allocate(coordsx(Nx,Nx),coordsy(Nx,Nx),forcing(Nx,Nx), &
	&		 streamfunc(Nx,Nx),friction(Nx,Nx),nonlinterm(Nx,Nx))

!!!!!!!!!!!
	ppi=3.141592653589793
	e=2.718281828459045

!for debug
	ltrace=.true.
	


!!!!!!!!!!!!! Make friction
if (fric==1) then
		friction=u(ivort,:,:)

	if (filterfric==1) then	!Apply Fourier filter if requested...
!		call fftfilter(friction,Nx,0,kfric)
	endif
endif
!Dump fricx/y and forcingx/y to file, to check out...
if(1==0)then
        if (mod(ITE-1,freq).eq.0) then
	if(printforcing.eq.1)then
	ret = gft_out_brief('fric',time, (/nx,nx/), 2, friction)
	endif
	endif
endif
        if (mod(ITE-1,freq).eq.0) then
	if(printforcing.eq.1)then
!	ret = gft_out_brief(trim(arg)//'f',time, (/nx,nx/), 2, forcing)
	endif
	endif

!Make edges match up
	forcing(:,Nx)=forcing(:,1)
	forcing(Nx,:)=forcing(1,:)
	friction(:,Nx)=friction(:,1)
	friction(Nx,:)=friction(1,:)

!!!!!!!!!!!
	source = 0.0
!	write(*,*) lapvort

	nonlinterm=-dxU(ivort,:,:)*dypsi(ivort,:,:) +&
		  & dyU(ivort,:,:)*dxpsi(ivort,:,:)

	call fftfilter(nonlinterm,Nx,msk,FFTHandle,lengths)

!if(round==1)then
	source(ivort,:,:) = nonlinterm &
		& + nu*lapvort(ivort,:,:) &
		& - fricstr*friction
!elseif(round==2)then

!endif
	
	deallocate(coordsx,coordsy,forcing,streamfunc,friction,nonlinterm)
			
	end subroutine rhs

	end module M_RHS
