MODULE POISSON_SOLVER
	use MKL_DFTI
	IMPLICIT NONE

CONTAINS

	subroutine inverse_laplacian(vort,Nx,dx,kx,ky,FFTHandle,lengths)
	use m_pars, only : Lx
	implicit none

	!Incoming objects
	integer :: Nx !grid size
	real*8 :: dx !spatial step size
	real*8, dimension(:,:), intent(inout) :: vort !vorticity field
	real*8, dimension(:,:), intent(in) :: kx,ky
	integer, dimension(:), intent(in) :: lengths

	!constants
	real :: pi=3.14159265359

	!For do loops
	integer :: i,j,l

	!MKL FFT stuff
    integer, parameter :: WP = selected_real_kind(15,307)
    real(WP), allocatable :: xx_real(:,:)
    real(WP), allocatable :: xx_real2(:,:)!, kk(:,:)
    type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
    integer :: sStatus

	!Working arrays
	real*8, dimension(:,:), allocatable :: kk

    allocate(xx_real(lengths(1),lengths(2)))
    allocate(xx_real2(lengths(1),lengths(2)))
	allocate(kk(lengths(1),lengths(2)))
	
	l=lengths(1)

	!Initialize Cr with vort array, redundant edges trimmed off
	xx_real=vort(1:l,1:l)

	!Transform vort to Fourier space
    sStatus = DftiComputeForward(FFTHandle, xx_real(:,1), xx_real2(:,1) )

	!Solve for stream function in Fourier space
!	xx_real2(2:l,2:l) = 0.5*dx**2.*xx_real2(2:l,2:l)&
!					&/(cos(2*pi*kx(2:l,2:l)/l) + cos(2*pi*ky(2:l,2:l)/l) - 2.)
!	xx_real2(2:l,1)   = 0.5*dx**2.*xx_real2(2:l,1)&
!					&/(cos(2*pi*kx(2:l,1)/l) + cos(2*pi*ky(2:l,1)/l) - 2.)
!	xx_real2(1,2:l)   = 0.5*dx**2.*xx_real2(1,2:l)&
!					&/(cos(2*pi*kx(1,2:l)/l) + cos(2*pi*ky(1,2:l)/l) - 2.)
!	xx_real2(1,1)     = 0.0

	kk=(2.*pi/Lx)**2.*(kx**2.+ky**2.)
	kk(1,1)=1. !Omit 0 frequency term
	xx_real2   = xx_real2/kk
	xx_real2(1,1)   = 0.0	!Kill constant term

	!Transform stream function back to real space
    sStatus = DftiComputeBackward(FFTHandle, xx_real2(:,1), xx_real(:,1) )

	!Dump the result into vort array, not including far edges
	vort(1:Nx-1,1:Nx-1)=xx_real

	!Make edges match up
	vort(:,Nx)=vort(:,1)
	vort(Nx,:)=vort(1,:)
	deallocate(xx_real,xx_real2,kk)
	!done
	end subroutine inverse_laplacian

END MODULE POISSON_SOLVER
