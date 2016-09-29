MODULE decompose
	IMPLICIT NONE

CONTAINS
	subroutine dcmpose(v,vavg,vfluc,Nx,width)
	IMPLICIT NONE
	INTEGER, INTENT(IN) :: Nx, width !width must be odd
	INTEGER :: i,j,m,n,i2,j2
	REAL, DIMENSION(:,:), INTENT(INOUT) :: vavg, vfluc
	REAL, DIMENSION(:,:), INTENT(IN) :: v
	REAL, DIMENSION(Nx,Nx) :: blotch
!	ALLOCATE(blotch(-width/2:width/2,-width/2:width/2))
	vavg=0.
	do m=1,Nx
		do n=1,Nx
!			vavg(i,j)=smooth(i,j,width,Nx)
			do i=0,width
			i2=wrappedi(m+i-width/2,Nx)
				do j=0,width
					j2=wrappedj(n+j-width/2,Nx)
					vavg(m,n)=vavg(m,n)+v(i2,j2)
				enddo
			enddo
		enddo
	enddo
	vfluc=v-vavg

	print *, 'made it thru.'
!	deallocate(blotch)
	CONTAINS

	INTEGER FUNCTION wrappedi(i,N)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: i	!coords of 'field point'.
		INTEGER :: N
		REAL :: newi
		N=N-1
		if (i<1) then
			newi=N+i
			wrappedi=newi
		elseif (i>N) then
			newi=i-N
			wrappedi=newi
		else
			wrappedi=i
		endif

	END FUNCTION wrappedi

	INTEGER FUNCTION wrappedj(j,N)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: j	!coords of 'field point'.
		INTEGER :: N
		REAL :: newj
		N=N-1
		if (j<1) then
			newj=N+j
			wrappedj=newj
		elseif (j>N) then
			newj=j-N
			wrappedj=newj
		else
			wrappedj=j
		endif
	END FUNCTION wrappedj

	REAL FUNCTION smooth(x,y,width,Nx)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: Nx
!		REAL, INTENT(IN), DIMENSION(Nx,Nx) :: v
!		REAL, DIMENSION(:,:), ALLOCATABLE :: blotch
		INTEGER, INTENT(IN) :: width
		REAL :: temp
		INTEGER :: x,y,i,j,i2,j2
!		ALLOCATE(blotch(-width/2:width/2,-width/2:width/2))
		temp=0.
		do i=0,width
			i2=wrappedi(int(x)+i-width/2,Nx)
			do j=0,width
				j2=wrappedj(int(y)+j-width/2,Nx)
				temp=temp+v(i2,j2)
			enddo
		enddo
		smooth=temp/float(width**2)
!		deallocate(blotch)
	END FUNCTION smooth
	end subroutine dcmpose
END MODULE decompose
