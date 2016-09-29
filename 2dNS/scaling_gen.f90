    module scaling_gen
    implicit none
!This uses the PACK packing convention, see the figure here:
!https://software.intel.com/en-us/node/470850#8217B6CE-36F4-463A-8E79-B88AF67FDA27

    contains


    subroutine kxky(kx,ky,N)
    implicit none

    real*8, dimension(:,:), intent(inout) :: kx,ky
    integer, intent(in) :: N
	integer :: M,K,L,i,j

	M=N

!!!!!! N even
	if (mod(N,2)==0) then
		K=M/2
		L=N/2
!!!kx
		!First row (red)
		kx(1,:)=0.

		!First and last columns
		do i=2,2*K-2,2
			!First column (green)
			kx(i,1)  =i/2.
			kx(i+1,1)=i/2.
			!Last column (blue)
			kx(i,2*L)  =i/2.
			kx(i+1,2*L)=i/2.
		enddo
		kx(2*K,1)  =K !Last element green column
		kx(2*K,2*L)=K !Last element blue column

		!Centre block (yellow)
		do i=2,2*K
			kx(i,2:2*L-1)=i-1.
		enddo

!!!ky
		!First column (green) and R_{0,0}
		ky(:,1)  =0.

		!Last column (blue) and R_{0,L}
		ky(:,2*L)=L

		!Everywhere else
		do i=2,2*L-2,2
			ky(:,i)  =i/2.
			ky(:,i+1)=i/2.
		enddo

!!!!!! N odd
	elseif (mod(N,2)==1) then
		K=(M-1)/2
		L=(N-1)/2

!!!kx
		!First row (red)
		kx(1,:)=0.

		!First column (green)
		do i=2,2*K,2
			kx(i,1)  =i/2.
			kx(i+1,1)=i/2.
		enddo

		!Yellow block
		do i=2,2*K+1
			kx(i,2:2*L+1)=i-1.
		enddo

!!!ky
		!First column (green) and R_{0,0}
		ky(:,1)=0.

		!Everywhere else
		do i=2,2*L,2
			ky(1:2*K+1,i)  =i/2.
			ky(1:2*K+1,i+1)=i/2.
		enddo

	endif

!Turn wavenumbers>Nyquist into i-N values instead
	do i=1,2*L
		do j=1,2*K
			if (kx(i,j)>M/2) then
				kx(i,j)=kx(i,j)-M
			endif
			if (ky(i,j)>N/2) then
				ky(i,i)=ky(i,j)-N
			endif
		enddo
	enddo

!Take transpose for fortran
	kx=TRANSPOSE(kx)
	ky=TRANSPOSE(ky)

	end subroutine kxky

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine real_imag_masks(kreal,kimag,N)
	implicit none

	logical, dimension(:,:), intent(inout) :: kreal, kimag
	integer, intent(in) :: N
	integer :: M,K,L,i,j

!!!!!! N even
        if (mod(N,2)==0) then
                K=M/2
                L=N/2
!!!kreal
		kreal=.False. !initialize

		kreal(1,1)     =.True. !Top left corner
		kreal(1,2*L)   =.True. !Top right corner
		kreal(2*K,1)   =.True. !Bottom left corner
		kreal(2*K,2*L) =.True. !Bottom right corner

		do i=2,2*L-2,2
			!All rows
			kreal(:,i)   =.True.
			!Left column (green)
			kreal(i,1)   =.True.
			!Right column (blue)
			kreal(i,2*L) =.True.
		enddo			

	endif
!!!!!! N odd
	if (mod(N,2)==1) then
		K=(M-1)/2
		L=(N-1)/2

		kreal=.False. !initialize

		kreal(1,1) =.True. !Top left corner
	
		do i=2,2*L,2
			!All but left column
			kreal(:,i) =.True.
			!Left column (green)
			kreal(i,1) =.True.
		enddo
	endif

!!!kimag
		kimag=.not.kreal

!Take transpose for fortran
        kreal=TRANSPOSE(kreal)
        kimag=TRANSPOSE(kimag)

        end subroutine real_imag_masks



	

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine rescale(scaling,F0)
	use MKL_DFTI
	use m_pars, only: N
	implicit none

	real*8, dimension(:,:), intent(inout) :: scaling
	real*8, intent(in) :: F0

    !MKL DFT stuff
	integer, parameter :: WP = selected_real_kind(15,307)
    real(WP), allocatable :: x_real(:,:)
	real(WP), allocatable :: x_real2(:,:)
	integer :: cstrides(3), rstrides(3)
    type(DFTI_DESCRIPTOR), POINTER :: My_Desc_Handle
    integer :: sStatus, lengths(2)
!	integer :: ret, gft_out_brief

	My_Desc_Handle => null()

    lengths(1)=N-1
    lengths(2)=N-1

	allocate(x_real(lengths(1),lengths(2)))
	allocate(x_real2(lengths(1),lengths(2)))

	x_real2=2.*scaling**2.

    sStatus = DftiCreateDescriptor(My_Desc_Handle, DFTI_DOUBLE,&
            & DFTI_REAL, 2, lengths)
	sStatus = DftiSetValue(My_Desc_Handle, DFTI_PLACEMENT,&
			& DFTI_NOT_INPLACE)
    sStatus = DftiSetValue(My_Desc_Handle, DFTI_CONJUGATE_EVEN_STORAGE,&
            & DFTI_COMPLEX_REAL)
    sStatus = DftiSetValue(My_Desc_Handle, DFTI_PACKED_FORMAT,&
            & DFTI_PACK_FORMAT)

	rstrides = [0, 1, lengths(1)]
	cstrides = [0, 1, lengths(2)]

	sStatus = DftiSetValue(My_Desc_Handle, DFTI_INPUT_STRIDES, rstrides)
	sStatus = DftiSetValue(My_Desc_Handle, DFTI_OUTPUT_STRIDES, cstrides)

    sStatus = DftiCommitDescriptor(My_Desc_Handle)

	sStatus = DftiComputeBackward(My_Desc_Handle, x_real2(:,1), x_real(:,1) )

!Recale 'scaling' to have desired (1,1) component in real space
    scaling=scaling*sqrt(F0/x_real(1,1))

!check
	x_real2=2.*scaling**2.
	sStatus = DftiComputeBackward( My_Desc_Handle, x_real2(:,1), x_real(:,1) )
    write(*,*) x_real(1,1)
!	ret = gft_out_brief('scaling',1.0d0,(/800,800/),2,x_real)
	sStatus = DftiFreeDescriptor(My_Desc_Handle)
	deallocate(x_real,x_real2)
	end subroutine rescale
	end module
