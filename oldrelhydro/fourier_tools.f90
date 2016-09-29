MODULE FOURIER_TOOLS
	IMPLICIT NONE

CONTAINS

	subroutine fftfilter(R,Nx,k1,k2)
	use fourier_repackage

	implicit none

	!Incoming objects
	integer :: Nx,k1,k2
	real*8, dimension(:,:) :: R

	!For CFFT2I or RFFT2I
	integer ( kind = 4 ) :: l
	integer ( kind = 4 ) :: m
	integer ( kind = 4 ) :: lensav
	integer ( kind = 4 ) :: ier
	real ( kind = 4 ), dimension(:), allocatable :: wsave
	
	!For CFFT2F/B
	integer ( kind = 4 ) :: LDIM
	integer ( kind = 4 ) :: LENWRK
	real ( kind = 4 ), dimension(:), allocatable :: WORK
	complex (kind=4), dimension(:,:), allocatable :: Cr

	!For do loops
	integer :: i,j
	real(kind=4) :: kmag

	!Coordinate arrays
	real ( kind = 4 ), dimension(:,:), allocatable :: coordx,coordy

	l=Nx
	m=Nx
	lensav=2*(L+M+INT(LOG(REAL(L)))+INT(LOG(REAL(M)))+8)
	LDIM=2*(l/2+1)
	LENWRK=2*LDIM*M   !LDIM*M if doing real array
	allocate(wsave(lensav),Cr(LDIM,M),WORK(LENWRK),coordx(L,M),coordy(L,M))

	!Define coordinates.
	do i=1,L
		coordy(i,:) = float(i-1)
		coordx(:,i) = float(i-1)
	enddo
	
	coordx=coordx-L/2	!Centre coords on the origin (this works for even and odd L..
	coordy=coordy-L/2	!Draw a picture if you don't believe it).

	!Initialize Cr
	Cr=cmplx(0.,0.)

	!Put data into Cr
	do j=1,Nx
		do i=1,Nx
			Cr(i,j)=cmplx(R(i,j),0.)
		enddo
	enddo

	!Transform to Fourier space
	call cfft2i ( l, m, wsave, lensav, ier )
	call cfft2f ( ldim, l, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Repackage into centre-origin form
	call repackagef(Cr,L,M)

	!Bin the appropriate entries to zero.
	do j=1,M
		do i=1,L
			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 ) 
			if ((kmag.lt.k1).or.(kmag.ge.k2)) then
				Cr(i,j)=cmplx(0.,0.)
			endif
		enddo
	enddo

	!Repackage into corner-origin form
	call repackageb(Cr,L,M)

	!Impose reality of the data (should already be real, but just in case)
	do j=1,M
		do i=1,L/2+1
			if((i==1).and.(j==1))then
				cycle
			elseif(i==1)then
				Cr(i,j) = conjg( Cr(i,M-(j-2)) )
			elseif(j==1)then
				Cr(i,j) = conjg( Cr(L-(i-2),j) )
			else
				Cr(i,j) = conjg( Cr(L-(i-2),M-(j-2)) )
			endif
		enddo
	enddo

	!Now transform back to real space
	call cfft2i ( l, m, wsave, lensav, ier )
	call cfft2b ( ldim, l, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=real(Cr(1:Nx-1,1:Nx-1))

	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(wsave,Cr,WORK,coordx,coordy)
	!done
	end subroutine fftfilter







	subroutine forcingfield(R,Nx,k1,k2,F0)
	use fourier_repackage

	implicit none

	!Incoming objects
	integer :: Nx,k1,k2
	real*8, dimension(:,:) :: R
	real*8 :: F0

	!For CFFT2I or RFFT2I
	integer :: eNx
	integer ( kind = 4 ) :: l
	integer ( kind = 4 ) :: m
	integer ( kind = 4 ) :: lensav
	integer ( kind = 4 ) :: ier
	real ( kind = 4 ), dimension(:), allocatable :: wsave
	
	!For CFFT2F/B
	integer ( kind = 4 ) :: LDIM
	integer ( kind = 4 ) :: LENWRK
	real ( kind = 4 ), dimension(:), allocatable :: WORK
	complex (kind=4), dimension(:,:), allocatable :: Cr

	!For do loops
	integer :: i,j
	real*4, dimension(:,:), allocatable :: kmag
	logical, dimension(:,:), allocatable :: maskr

	!Coordinate arrays
	real ( kind = 4 ), dimension(:,:), allocatable :: coordx,coordy

	!For random number generator
	integer :: idum

	l=Nx
	m=Nx
	lensav=2*(L+M+INT(LOG(REAL(L)))+INT(LOG(REAL(M)))+8)
	LDIM=2*(l/2+1)
	LENWRK=2*LDIM*M   !LDIM*M if doing real array
	allocate(wsave(lensav),Cr(LDIM,M),WORK(LENWRK),coordx(L,M),coordy(L,M),kmag(L,M),maskr(L,M))

	!Define coordinates.
	do i=1,L
		coordy(i,:) = float(i-1)
		coordx(:,i) = float(i-1)
	enddo
	
	coordx=coordx-L/2	!Centre coords on the origin (this works for even and odd L..
	coordy=coordy-L/2	!Draw a picture if you don't believe it).

	!Initialize Cr
	Cr=cmplx(0.,0.)

	!Fill Cr with zero-mean unit-variance (complex) Gaussian random numbers.
	do j=1,Nx
		do i=1,Nx
			Cr(i,j)=cmplx(sqrt(F0)*gasdev(idum),sqrt(F0)*gasdev(idum))
		enddo
	enddo

	!Bin the appropriate entries to zero.
	kmag=coordx**2 + coordy**2
	maskr=(kmag.ge.k1**2).and.(kmag.lt.k2**2)
!	do j=1,M
!		do i=1,L
!			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 ) 
!			if ((kmag(i,j).lt.k1).or.(kmag(i,j).ge.k2)) then
			where(.not.maskr) Cr=cmplx(0.,0.)
!			endif
!		enddo
!	enddo

	!Repackage into corner-origin form
	call repackageb(Cr,L,M)

	!Impose reality of the data (should already be real, but just in case)
	do j=1,M
		do i=1,L/2+1
			if((i==1).and.(j==1))then
				cycle
			elseif(i==1)then
				Cr(i,j) = conjg( Cr(i,M-(j-2)) )
			elseif(j==1)then
				Cr(i,j) = conjg( Cr(L-(i-2),j) )
			else
				Cr(i,j) = conjg( Cr(L-(i-2),M-(j-2)) )
			endif
		enddo
	enddo

	!Now transform back to real space
	call cfft2i ( l, m, wsave, lensav, ier )
	call cfft2b ( ldim, l, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=real(Cr(1:Nx-1,1:Nx-1))

	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(wsave,Cr,WORK,coordx,coordy)
	!done
	end subroutine forcingfield






	subroutine forcingfield2(R,Nx,lf,F0) !Forcing with Gaussian real-space corr function 
	use fourier_repackage

	implicit none

	!Incoming objects
	integer :: Nx, lf !lf is real-space correlation length in grid units.
	real*8 :: F0 !Energy injection rate.
	real*8, dimension(:,:) :: R

	!constants
	real :: pi=3.14159265359

	!For CFFT2I or RFFT2I
	integer :: eNx
	integer ( kind = 4 ) :: l
	integer ( kind = 4 ) :: m
	integer ( kind = 4 ) :: lensav
	integer ( kind = 4 ) :: ier
	real ( kind = 4 ), dimension(:), allocatable :: wsave
	
	!For CFFT2F/B
	integer ( kind = 4 ) :: LDIM
	integer ( kind = 4 ) :: LENWRK
	real ( kind = 4 ), dimension(:), allocatable :: WORK
	complex (kind=4), dimension(:,:), allocatable :: Cr

	!For do loops
	integer :: i,j
	real(kind=4) :: kmag
	real :: scaling,lf2

	!Coordinate arrays
	real ( kind = 4 ), dimension(:,:), allocatable :: coordx,coordy

	!For random number generator
	integer :: idum

	l=Nx-1
	m=Nx-1
	lensav=2*(L+M+INT(LOG(REAL(L)))+INT(LOG(REAL(M)))+8)
	LDIM=l
	LENWRK=2*LDIM*M   !LDIM*M if doing real array
	allocate(wsave(lensav),Cr(LDIM,M),WORK(LENWRK),coordx(L,M),coordy(L,M))

	!Define coordinates.
	do i=1,L
		coordy(i,:) = float(i-1)
		coordx(:,i) = float(i-1)
	enddo
	
	coordx=coordx-L/2	!Centre coords on the origin (this works for even and odd L..
	coordy=coordy-L/2	!Draw a picture if you don't believe it).

	!Initialize Cr
	Cr=cmplx(0.,0.)

	!Fill Cr with zero-mean unit-variance (complex) Gaussian random numbers.
	do j=1,m
		do i=1,ldim
			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 )
			kmag=2*pi*kmag/10.
			lf2=10.*lf/float(Nx)
			scaling=sqrt(F0*lf2**2*exp(-lf2**2*kmag**2/2.))
!			write(*,*) -2.*pi**2*lf**2*kmag**2
			Cr(i,j)=cmplx(scaling*gasdev(idum),scaling*gasdev(idum))
		enddo
	enddo

	!Bin the appropriate entries to zero.
!	do j=1,M
!		do i=1,L
!			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 ) 
!			Cr(i,j)=Cr(i,j)*sqrt(F0*pi*0.5*exp(-2.*pi**2*lf**2*kmag**2))
!		enddo
!	enddo

	!Repackage into corner-origin form
	call repackageb(Cr,Ldim,M)

	!Impose reality of the data (should already be real, but just in case)
	do j=1,M
		do i=1,Ldim/2+1
			if((i==1).and.(j==1))then
				cycle
			elseif(i==1)then
				Cr(i,j) = conjg( Cr(i,M-(j-2)) )
			elseif(j==1)then
				Cr(i,j) = conjg( Cr(L-(i-2),j) )
			else
				Cr(i,j) = conjg( Cr(L-(i-2),M-(j-2)) )
			endif
		enddo
	enddo

	!Now transform back to real space
	call cfft2i ( l, m, wsave, lensav, ier )
	call cfft2b ( ldim, l, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=real(Cr(1:Nx-1,1:Nx-1))

	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(wsave,Cr,WORK,coordx,coordy)
	!done
	end subroutine forcingfield2






	subroutine forcingfield3(R,Nx,w,s,F0) !Forcing with Gaussian shifted profile in Fourier space,
					      !and -> 0 quickly in real space.
	use fourier_repackage

	implicit none

	!Incoming objects
	integer :: Nx, w, s !lf is real-space correlation length in grid units.
	real*8 :: F0 !Energy injection rate.
	real*8, dimension(:,:) :: R

	!constants
	real :: pi=3.14159265359

	!For CFFT2I or RFFT2I
	integer :: eNx
	integer ( kind = 4 ) :: l
	integer ( kind = 4 ) :: m
	integer ( kind = 4 ) :: lensav
	integer ( kind = 4 ) :: ier
	real ( kind = 4 ), dimension(:), allocatable :: wsave
	
	!For CFFT2F/B
	integer ( kind = 4 ) :: LDIM
	integer ( kind = 4 ) :: LENWRK
	real ( kind = 4 ), dimension(:), allocatable :: WORK
	complex (kind=4), dimension(:,:), allocatable :: Cr

	!For do loops
	integer :: i,j
	real(kind=4) :: kmag
	real :: scaling,w2,s2

	!Coordinate arrays
	real ( kind = 4 ), dimension(:,:), allocatable :: coordx,coordy

	!For random number generator
	integer :: idum

	l=Nx
	m=Nx
	lensav=2*(L+M+INT(LOG(REAL(L)))+INT(LOG(REAL(M)))+8)
	LDIM=2*(l/2+1)
	LENWRK=2*LDIM*M   !LDIM*M if doing real array
	allocate(wsave(lensav),Cr(LDIM,M),WORK(LENWRK),coordx(L,M),coordy(L,M))

	!Define coordinates.
	do i=1,L
		coordy(i,:) = float(i-1)
		coordx(:,i) = float(i-1)
	enddo
	
	coordx=coordx-L/2	!Centre coords on the origin (this works for even and odd L..
	coordy=coordy-L/2	!Draw a picture if you don't believe it).

	!Initialize Cr
	Cr=cmplx(0.,0.)

	!Fill Cr with zero-mean unit-variance (complex) Gaussian random numbers, then scale
	!them to have the variance dictated by the 2pt correlation function.
	do j=1,Nx
		do i=1,Nx
			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 )
			kmag=2*pi*kmag/10.
			if(kmag==0.)then
				scaling=0.
			else
				w2=10.*w/float(Nx)	!In real space.
				s2=2*pi*s/10.		!In fourier space.
				scaling=sqrt(F0*w2**2*exp(-w2**2*(kmag-s2)**2/2.))
			endif

			if(kmag==0.)then
				Cr(i,j)=cmplx(scaling*gasdev(idum),scaling*gasdev(idum))
			else
				Cr(i,j)=cmplx(scaling*gasdev(idum),scaling*gasdev(idum))!/kmag**2.
			endif
		enddo
	enddo

	!Bin the appropriate entries to zero.
!	do j=1,M
!		do i=1,L
!			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 ) 
!			Cr(i,j)=Cr(i,j)*sqrt(F0*pi*0.5*exp(-2.*pi**2*lf**2*kmag**2))
!		enddo
!	enddo

	!Repackage into corner-origin form
	call repackageb(Cr,L,M)

	!Impose reality of the data (should already be real, but just in case)
	do j=1,M
		do i=1,L/2+1
			if((i==1).and.(j==1))then
				cycle
			elseif(i==1)then
				Cr(i,j) = conjg( Cr(i,M-(j-2)) )
			elseif(j==1)then
				Cr(i,j) = conjg( Cr(L-(i-2),j) )
			else
				Cr(i,j) = conjg( Cr(L-(i-2),M-(j-2)) )
			endif
		enddo
	enddo

	!Now transform back to real space
	call cfft2i ( l, m, wsave, lensav, ier )
	call cfft2b ( ldim, l, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=real(Cr(1:Nx-1,1:Nx-1))

	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(wsave,Cr,WORK,coordx,coordy)
	!done
	end subroutine forcingfield3




	subroutine forcingfield4(R,Nx,w,s,F0) !Forcing with shifed raised cosine profile in Fourier 						      !space, and -> 0 quickly in real space.
	use fourier_repackage

	implicit none

	!Incoming objects
	integer :: Nx, w, s !lf is real-space correlation length in grid units.
	real*8 :: F0 !Energy injection rate.
	real*8, dimension(:,:) :: R

	!constants
	real :: pi=3.14159265359

	!For CFFT2I or RFFT2I
	integer :: eNx
	integer ( kind = 4 ) :: l
	integer ( kind = 4 ) :: m
	integer ( kind = 4 ) :: lensav
	integer ( kind = 4 ) :: ier
	real ( kind = 4 ), dimension(:), allocatable :: wsave
	
	!For CFFT2F/B
	integer ( kind = 4 ) :: LDIM
	integer ( kind = 4 ) :: LENWRK
	real ( kind = 4 ), dimension(:), allocatable :: WORK
	complex (kind=4), dimension(:,:), allocatable :: Cr

	!For do loops
	integer :: i,j
	real(kind=4) :: kmag
	real :: scaling,w2,s2

	!Coordinate arrays
	real ( kind = 4 ), dimension(:,:), allocatable :: coordx,coordy,kmagsq,scal
    logical, dimension(:,:), allocatable :: maskr

	!For random number generator
	integer :: idum

	l=Nx
	m=Nx
	lensav=2*(L+M+INT(LOG(REAL(L)))+INT(LOG(REAL(M)))+8)
	LDIM=2*(l/2+1)
	LENWRK=2*LDIM*M   !LDIM*M if doing real array
	allocate(wsave(lensav),Cr(LDIM,M),WORK(LENWRK),coordx(L,M),coordy(L,M),kmagsq(L,M),scal(L,M),maskr(L,M))

	!Define coordinates.
	do i=1,L
		coordy(i,:) = float(i-1)
		coordx(:,i) = float(i-1)
	enddo
	
	coordx=coordx-L/2	!Centre coords on the origin (this works for even and odd L..
	coordy=coordy-L/2	!Draw a picture if you don't believe it).

	!Initialize Cr
	Cr=cmplx(0.,0.)

	!Fill Cr with zero-mean unit-variance (complex) Gaussian random numbers, then scale
	!them to have the variance dictated by the 2pt correlation function.
	do j=1,Nx
		do i=1,Nx
            Cr(i,j)=cmplx(gasdev(idum),gasdev(idum))
        enddo
    enddo

	kmagsq=(2*pi/10.)**2*(coordx**2+coordy**2)
	scal=0.	!initialize
	w2=(2*pi*w/10.)
	s2=(2*pi*s/10.)
	maskr=(kmagsq.ge.(s2-w2)**2.).and.(kmagsq.le.(s2+w2)**2.)
	where( maskr ) scal=sqrt(10*0.5*F0*(1/w2)**2*(1+cos((sqrt(kmagsq)-s2)*pi/w2)))

    Cr(1:Nx,1:Nx)=scal*Cr(1:Nx,1:Nx)

	!Repackage into corner-origin form
	call repackageb(Cr,L,M)

	!Impose reality of the data (should already be real, but just in case)
	do j=1,M
		do i=1,L/2+1
			if((i==1).and.(j==1))then
				cycle
			elseif(i==1)then
				Cr(i,j) = conjg( Cr(i,M-(j-2)) )
			elseif(j==1)then
				Cr(i,j) = conjg( Cr(L-(i-2),j) )
			else
				Cr(i,j) = conjg( Cr(L-(i-2),M-(j-2)) )
			endif
		enddo
	enddo

	!Now transform back to real space
	call cfft2i ( l, m, wsave, lensav, ier )
	call cfft2b ( ldim, l, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=real(Cr(1:Nx-1,1:Nx-1))

	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(wsave,Cr,WORK,coordx,coordy)
	!done
	end subroutine forcingfield4








	!This one works well.
	subroutine forcingfield5(R,Nx,lf,F0) !Forcing with Gaussian/k^2 Fourier spectrum 
	use fourier_repackage

	implicit none

	!Incoming objects
	integer :: Nx, lf !lf is real-space correlation length in grid units.
	real*8 :: F0 !Energy injection rate.
	real*8, dimension(:,:) :: R

	!constants
	real :: pi=3.14159265359

	!For CFFT2I or RFFT2I
	integer :: eNx
	integer ( kind = 4 ) :: l
	integer ( kind = 4 ) :: m
	integer ( kind = 4 ) :: lensav
	integer ( kind = 4 ) :: ier
	real ( kind = 4 ), dimension(:), allocatable :: wsave
	
	!For CFFT2F/B
	integer ( kind = 4 ) :: LDIM
	integer ( kind = 4 ) :: LENWRK
	real ( kind = 4 ), dimension(:), allocatable :: WORK
	complex (kind=4), dimension(:,:), allocatable :: Cr

	!For do loops
	integer :: i,j
	real(kind=4) :: kmag
	real :: scaling,lf2

	!Coordinate arrays
	real ( kind = 4 ), dimension(:,:), allocatable :: coordx,coordy

	!For random number generator
	integer :: idum

	l=Nx
	m=Nx-1
	lensav=2*(M+M+INT(LOG(REAL(M)))+INT(LOG(REAL(M)))+8)
	LDIM=l-1 !2*(l/2+1)
	LENWRK=2*LDIM*M   !LDIM*M if doing real array
	allocate(wsave(lensav),Cr(LDIM,M),WORK(LENWRK),coordx(Ldim,M),coordy(Ldim,M))

	!Define coordinates.
	do i=1,Ldim
		coordy(i,:) = float(i-1)
		coordx(:,i) = float(i-1)
	enddo
	
	coordx=coordx-Ldim/2	!Centre coords on the origin (this works for even and odd L..
	coordy=coordy-Ldim/2	!Draw a picture if you don't believe it - I think only works for
				! even L... April 22 2015).

	!Initialize Cr
	Cr=cmplx(0.,0.)

	!Fill Cr with zero-mean unit-variance (complex) Gaussian random numbers.
	do j=1,m
		do i=1,ldim
			kmag=sqrt( coordx(i,j)**2 + coordy(i,j)**2 )
			kmag=2*pi*kmag/10.
			lf2=10.*lf/float(Nx)
			if(kmag==0)then
			scaling=0.
			else
			scaling=sqrt(F0*lf2**2*exp(-lf2**2*kmag**2/2.))/kmag**2.
			endif
			Cr(i,j)=cmplx(scaling*gasdev(idum),scaling*gasdev(idum))
		enddo
	enddo

	!Repackage into corner-origin form
	call repackageb(Cr,Ldim,M)

	!Impose reality of the data. This assumes even dimensional array.
	do j=1,M
		do i=1,Ldim/2+1
			if( ((i==1).and.(j==1)).or.((i==1).and.(j==M-M/2+1)).or. &
				& ((i==Ldim-Ldim/2+1).and.(j==1)).or. &
				& ((i==Ldim-Ldim/2+1).and.(j==M-M/2+1)) )then
				Cr(i,j) = cmplx(real(Cr(i,j)),0.)
			elseif(i==1)then
				Cr(i,j) = conjg( Cr(i,M-(j-2)) )
			elseif(j==1)then
				Cr(i,j) = conjg( Cr(Ldim-(i-2),j) )
			else
				Cr(i,j) = conjg( Cr(Ldim-(i-2),M-(j-2)) )
			endif
		enddo
	enddo

	!Now transform back to real space
	call cfft2i ( ldim, m, wsave, lensav, ier )
	call cfft2b ( ldim, ldim, m, Cr, wsave, lensav, work, lenwrk, ier )

	!Dump back into R, leaving edge empty
	R(1:Nx-1,1:Nx-1)=real(Cr(1:Nx-1,1:Nx-1))
	!Make edges match up
	R(:,Nx)=R(:,1)
	R(Nx,:)=R(1,:)
	deallocate(wsave,Cr,WORK,coordx,coordy)
	!done
	end subroutine forcingfield5










	REAL FUNCTION gasdev(idum)		!Gives standard gaussian distributed random numbers
		integer :: iset,idum
		real :: fac,gset,rsq,v1,v2,ran1
		save iset,gset
		data iset/0/
		if (iset.eq.0) then
			v1=1.
			v2=1.
			rsq=v1**2+v2**2
			do while (rsq.ge.1..or.rsq.eq.0.)
				call random_number(v1)
				call random_number(v2)
				v1=2.*v1-1.
				v2=2.*v2-1.
				rsq=v1**2+v2**2
			enddo
			fac=sqrt(-2.*log(rsq)/rsq)
			gset=v1*fac
			gasdev=v2*fac
			iset=1
		else
			gasdev=gset
			iset=0
		endif
		return
	END FUNCTION gasdev

END MODULE FOURIER_TOOLS
