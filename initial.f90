	MODULE M_INITIAL
	implicit none
	
	CONTAINS
	
	subroutine initial(Uold,Pold,dx,dy,x,y,Nv,N)
	use m_pars, only : irho, ivx,ivy,resume,resume_integer,weos
	real*8, dimension(:,:,:) :: Uold, Pold
	real*8, dimension(:,:) :: x,y
	real*8 :: dx, dy
	integer :: Nv, N, i, ret, gft_read_brief, tslice, j
	
	real*8 :: ppi
	real*8, dimension(:,:), allocatable :: WW, WW2
	real*8, dimension(:,:,:), allocatable :: perturbation
	integer :: STAT, ierr
!	call init_random_seed()
!	call random_number(perturbation)	

!	perturbation(1,:,1)=perturbation(1,:,N)
!	perturbation(1,1,:)=perturbation(1,N,:)

	ppi = 3.141592653589793

	allocate(WW(N,N),STAT=ierr)
	allocate(WW2(N,N),STAT=ierr)
	allocate(perturbation(2,N,N))

	uold(:,:,:) = 0.0d0
	pold(:,:,:) = 0.0d0

	pold(irho,:,:) = 1.

if (resume==0) then
	do i=1,N
!	    if(y(1,i).gt.3.and.y(1,i).lt.7.) then
!	      pold(ivx,:,i) =  0.2*sin(20.*ppi/y(3,N)*y(1,i))!* &
!     &                                (y(1,i)-3.)**2*(y(1,i)-7.)**2/16.
!	    end if
!          if(x(i,1).gt.3.and.x(i,1).lt.7.) then
!              pold(ivy,i,:) =  0.5*sin(5.*ppi/x(N,3)*x(i,1))* &
!     &                                (x(i,1)-3.)**2*(x(i,1)-7.)**2/16.
!
!	    end if
	end do

!	pold(ivx,:,:)=pold(ivx,:,:) + 0.02*perturbation(1,:,:)
endif

!"Resume" block.
if (resume==1) then
	tslice=resume_integer
	ret = gft_read_brief('eps',tslice,pold(irho,:,:))
	ret = gft_read_brief('vx',tslice,pold(ivx,:,:))
	ret = gft_read_brief('vy',tslice,pold(ivy,:,:))

endif


!now, get the conservatives
	WW = sqrt(1.-pold(ivx,:,:)**2 - pold(ivy,:,:)**2)
	WW2 = (1.-pold(ivx,:,:)**2-pold(ivy,:,:)**2)

	uold(irho,:,:) = pold(irho,:,:)*(-1.0d0+(weos+1.)/WW2)
	uold(ivx,:,:) = (weos+1.)*pold(irho,:,:)*pold(ivx,:,:)/WW2
	uold(ivy,:,:) = (weos+1.)*pold(irho,:,:)*pold(ivy,:,:)/WW2

	deallocate(WW,WW2,perturbation)
	
	end subroutine initial

	END MODULE M_INITIAL
