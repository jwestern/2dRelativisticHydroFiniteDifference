PROGRAM randomforcing
	use omp_lib
	IMPLICIT NONE
	INTEGER, PARAMETER :: Nx=401
	REAL, DIMENSION(Nx,Nx) :: forcingx, coordsx, coordsy, coordsx_rot, coordsy_rot
	REAL :: x,y,k,k2
	INTEGER :: kx, ky, i, j
	REAL, PARAMETER :: ppi=3.14159265359
!if (1==0) then
!Define array coordinates
	do i=1,Nx
		coordsx(i,:)=i
	enddo
	do j=1,Nx
		coordsy(:,j)=j
	enddo
	coordsx_rot=coordsx*cos(.3*ppi/2.)-coordsy*sin(.3*ppi/2.)
	coordsy_rot=coordsx*sin(.3*ppi/2.)+coordsy*cos(.3*ppi/2.)

	k=2*ppi*30./10.
	forcingx = sin(k*coordsx_rot/40.)*sin(k*coordsy_rot/40.)

		 !((cos(k*coordsx/40.)*sin(k*coordsy/40.)-sin(k*coordsx/40.)*cos(k*coordsy/40.))**2 &
			!& + (-cos(k*coordsx/40.)*sin(k*coordsy/40.)-sin(k*coordsx/40.)*cos(k*coordsy/40.))**2-1.0)
!	k2=25.
!	forcingx=forcingx + ((cos(k2*coordsx/40.)*sin(k2*coordsy/40.)-sin(k2*coordsx/40.)*cos(k2*coordsy/40.))**2 &
!			& + (-cos(k2*coordsx/40.)*sin(k2*coordsy/40.)-sin(k2*coordsx/40.)*cos(k2*coordsy/40.))**2-1.0)

if (1==0) then
do i=1,Nx
	x=10.*float(i)/float(Nx)
	do j=1,Nx
		y=10.*float(j)/float(Nx)	!Turn this into a proper coordinate (L=10)
		forcingx(i,j)=sin(60.*y)*sin(60.*x)
	enddo
enddo
print *, forcingx(1:5,1)
endif
if (1==0) then
!Test to see if openmp is working.
	!$omp parallel
!		print *, 'Hello world from thread', omp_get_thread_num(), '!'
	!$omp end parallel

!Define array coordinates
	do i=1,Nx
		coordsx(i,:)=i
	enddo
	do j=1,Nx
		coordsy(:,j)=j
	enddo

!Anders type forcing
	forcingx=0.
!$omp parallel private(kx,ky) shared(forcingx)
!$omp do
	do kx=-Nx/2,Nx/2-1
		do ky=-Nx/2,Nx/2-1
		if (omp_get_thread_num()==0) then
			print *, 'kx=', kx, 'ky=', ky
		endif
forcingx=forcingx+G(kx,ky)*(cos(2.*ppi*kx*coordsx/400.)*cos(2.*ppi*ky*coordsy/400.) &
					&- sin(2.*ppi*kx*coordsx/400.)*sin(2.*ppi*ky*coordsy/400.))
		enddo
	enddo
!$omp end do
!$omp end parallel
endif
	OPEN(1,FILE='sinkxsinkyforcing.dat',STATUS='REPLACE',ACCESS='SEQUENTIAL',ACTION='WRITE')
	WRITE(1,*) forcingx

CONTAINS

	REAL FUNCTION G(kx,ky)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: kx, ky
		if ((sqrt(float(kx**2+ky**2))<100.).and.(sqrt(float(kx**2+ky**2))>80.)) then
			G=1.
		else
			G=0.
		endif
	END FUNCTION G

END PROGRAM randomforcing
