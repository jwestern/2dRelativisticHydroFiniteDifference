MODULE FOURIER_TOOLS
	use MKL_VSL_TYPE
	use MKL_VSL
	use MKL_DFTI
	IMPLICIT NONE

CONTAINS

!Forcing with Gaussian real-space corr function 

	subroutine forcingfield2(R,Nx,scaling,errcode,method,stream,numrands,avg,sigma,FFTHandle,lengths)

	implicit none

	!Incoming objects
	integer :: Nx, lf !lf is real-space correlation length in grid units.
	real*8, dimension(:,:), intent(inout) :: R
	real*8, dimension(:,:), intent(in) :: scaling
    !For VSL RNG
    real(kind=8), intent(in) :: avg, sigma ! parameters of normal distribution
    TYPE (VSL_STREAM_STATE), intent(inout) :: stream
    integer(kind=4), intent(inout) :: errcode
    integer, intent(in) :: method,numrands
	integer, dimension(:), intent(in) :: lengths

	real*8, dimension(:), allocatable :: buff

	!constants
	real :: pi=3.14159265359

	!For do loops
	integer :: i,j,cnt

	integer, parameter :: WP = selected_real_kind(15,307)
    real(WP), allocatable :: x_real(:,:)
    real(WP), allocatable :: x_real2(:,:)
    type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
    integer :: sStatus

    allocate(x_real(lengths(1),lengths(2)))
    allocate(x_real2(lengths(1),lengths(2)))
	allocate(buff(numrands))

	!Generate all random numbers in buffer
	errcode=vdrnggaussian( method, stream, numrands, buff, avg, sigma )

	!Fill array with random numbers
	cnt=1
	do j=1,lengths(1)
		do i=1,lengths(2)
			x_real2(i,j)=buff(cnt)
			cnt=cnt+1
		enddo
	enddo

	!Scale the gaussian distributions at each point in fourier space
	x_real2=x_real2*scaling

	!Transform to real space
    sStatus = DftiComputeBackward(FFTHandle, x_real2(:,1), x_real(:,1) )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=x_real(1:Nx-1,1:Nx-1)
	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(x_real,x_real2,buff)
	!done
	end subroutine forcingfield2


	subroutine fftfilter(R,Nx,msk,FFTHandle,lengths)
	
        implicit none

        !Incoming objects
        integer :: Nx, lf !lf is real-space correlation length in grid units.
        real*8, dimension(:,:), intent(inout) :: R
        real*8, dimension(:,:), intent(in) :: msk
	integer, dimension(:), intent(in) :: lengths	

        !constants
        real :: pi=3.14159265359

        integer, parameter :: WP = selected_real_kind(15,307)
	real(WP), allocatable :: x_real(:,:)
	real(WP), allocatable :: x_real2(:,:)
	type(DFTI_DESCRIPTOR), POINTER, intent(inout) :: FFTHandle
	integer :: sStatus

	allocate(x_real(lengths(1),lengths(2)))
	allocate(x_real2(lengths(1),lengths(2)))

	x_real=R(1:Nx-1,1:Nx-1)

	sStatus = DftiComputeForward(FFTHandle, x_real(:,1), x_real2(:,1) )

	x_real2=x_real2*msk
	
	sStatus = DftiComputeBackward(FFTHandle, x_real2(:,1), x_real(:,1) )

	R(1:Nx-1,1:Nx-1)=x_real
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)

	deallocate(x_real,x_real2)

	end subroutine fftfilter

END MODULE FOURIER_TOOLS
