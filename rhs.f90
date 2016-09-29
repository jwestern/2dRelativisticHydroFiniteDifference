 	MODULE M_RHS
	implicit none
	
	contains
	
	subroutine rhs(u,dxu,dyU,dxFxU,dyFyU,source,dx,Nx,pold,FFTHandle,lengths)

	use m_pars, only : c, s,k_r, k_e, gamma,irho,ivx,ivy,solenoidalforce,fric,relativisticfriction, &
			& fricstr, shiftrectforce, kforcel, kforceg, kfric, filterfric,lf,F0, &
			& gaussforce,weos
	use fourier_tools
	use M_derivs

	implicit none
	real*8, dimension(:,:,:), intent(in) :: u, dxu,dyu, dxFxU,dyFyU,pold
	real*8, dimension(:,:,:), intent(inout) :: source
	type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
	integer, dimension(:), intent(inout) :: lengths
	real*8 :: dx, ppi, e 
	integer :: i, j, Nx, ret, gft_out_brief
	logical :: ltrace
	real*8, dimension(:,:), allocatable :: fricx,fricy

	allocate(fricx(Nx,Nx),fricy(Nx,Nx))

!!!!!!!!!!!
	ppi=3.141592653589793
	e=2.718281828459045

!for debug
	ltrace=.true.
	
!!!!!!!!!!!!! Make friction
if (fric==1) then
	if (relativisticfriction==1) then
		fricx=fricstr*u(ivx,:,:)
		fricy=fricstr*u(ivy,:,:)
	elseif(relativisticfriction==0)then
		fricx=fricstr*pold(ivx,:,:)
		fricy=fricstr*pold(ivy,:,:)
	endif

!	if (filterfric==1) then	!Apply Fourier filter if requested...
!XXX: NEED TO BE UPDATED TO USE MODERN FFTFILTER SUBROUTINE
!		call fftfilter(fricx,Nx,0,kfric)
!		call fftfilter(fricy,Nx,0,kfric)
!	endif
endif

!Make edges match up, just in case
	fricx(:,Nx)=fricx(:,1)
	fricx(Nx,:)=fricx(1,:)
	fricy(:,Nx)=fricy(:,1)
	fricy(Nx,:)=fricy(1,:)



!!!!!!!!!!!
	source = 0.0
		
	source(irho,:,:) = -dxFxU(irho,:,:) - dyFyU(irho,:,:)

	source(ivx,:,:) = -dxFxU(ivx,:,:) - dyFyU(ivx,:,:) - fricx
	
	source(ivy,:,:) = -dxFxU(ivy,:,:) - dyFyU(ivy,:,:) - fricy

	deallocate(fricx,fricy)
			
	end subroutine rhs
	
	end module M_RHS
