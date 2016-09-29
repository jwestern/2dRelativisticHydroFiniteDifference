MODULE FOURIER_REPACKAGE
	IMPLICIT NONE

CONTAINS

	SUBROUTINE repackagef(array,L,M) !From convoluted form to normal form
		implicit none

		complex (kind=4), dimension(:,:), allocatable :: work
		complex (kind=4), dimension(:,:) :: array
		integer :: L, M
		
		allocate(work(L,M))!array(2*(L/2+1),M)
		
		if(L.ne.M) print*, "Warning: not a square array... May not work..."

		if(mod(L,2).ne.0) then		!Check to see if dim is odd.
			work(L/2+1:L,M/2+1:M) = array(1:L/2+1,1:M/2+1)	!Move UL block to LR.
			work(1:L/2,1:M/2) = array(L/2+2:L,M/2+2:M)	!Move LR block to UL.
			work(L/2+1:L,1:M/2) = array(1:L/2+1,M/2+2:M)	!Move UR block to LL.
			work(1:L/2,M/2+1:M) = array(L/2+2:L,1:M/2+1)	!Move LL block to UR.

		else	!If it is even, then...
			work(L/2+1:L,M/2+1:M) = array(1:L/2,1:M/2)		!Move UL block to LR.
			work(1:L/2,1:M/2) = array(L/2+1:L,M/2+1:M)		!Move LR block to UL.
			work(L/2+1:L,1:M/2) = array(1:L/2,M/2+1:M)		!Move UR block to LL.
			work(1:L/2,M/2+1:M) = array(L/2+1:L,1:M/2)		!Move LL block to UR.
		endif
		array=0.
		array(1:L,1:M)=work
		deallocate(work)
	END SUBROUTINE repackagef

	SUBROUTINE repackageb(array,L,M) !From normal form to convoluted form
		implicit none

		complex (kind=4), dimension(:,:), allocatable :: work
		complex (kind=4), dimension(:,:) :: array
		integer :: L, M
		
		allocate(work(L,M))!array(2*(L/2+1),M)

		if(L.ne.M) print*, "Warning: not a square array... May not work..."

		if(mod(L,2).ne.0) then		!Check to see if dim is odd.
			work(1:L/2+1,1:M/2+1) = array(L/2+1:L,M/2+1:M)	!Move LR block to UR.
			work(L/2+2:L,M/2+2:M) = array(1:L/2,1:M/2)	!Move UL block to LR.
			work(1:L/2+1,M/2+2:M) = array(L/2+1:L,1:M/2)	!Move LL block to UR.
			work(L/2+2:L,1:M/2+1) = array(1:L/2,M/2+1:M)	!Move UR block to LL.
		else
			work(1:L/2,1:M/2) = array(L/2+1:L,M/2+1:M)		!Move UL block to LR.
			work(L/2+1:L,M/2+1:M) = array(1:L/2,1:M/2)		!Move LR block to UL.
			work(1:L/2,M/2+1:M) = array(L/2+1:L,1:M/2)		!Move UR block to LL.
			work(L/2+1:L,1:M/2) = array(1:L/2,M/2+1:M)		!Move LL block to UR.
		endif
		array=0.
		array(1:L,1:M)=work
	deallocate(work)
	END SUBROUTINE repackageb
END MODULE FOURIER_REPACKAGE
