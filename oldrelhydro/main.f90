	program main
	use M_initial
	use M_evolve
	use M_pars
	use M_derivs
	
	implicit none
	
!U will be conservative vars
!P will be primitive vars
	real*8, dimension(:,:,:), allocatable :: Unew, Pnew	
	real*8, dimension(:,:,:), allocatable :: Uold,Pold	
	real*8, dimension(:,:), allocatable :: x, y
 	real*8, dimension(:,:), allocatable :: vx_y, vy_x!, fx, fy
	real*8 :: dx,dy, dt, time
	integer :: ITE, i,j, Nv, Nsteps
    integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
	logical :: ltrace
	integer :: STAT,ierr
	character(len=32) :: arg

!Read in process number from command line
	CALL getarg(1,arg)

!set for debug or not
	ltrace = .false.

!read pars
        call readpars


!set the max number of eqns for now
	Nv = 3

!allocate needed memmory
!Unew will have the updated values of the fields
!uold will have the old values of the fields
!x, y are coords 
! var 1 --> D or e (cons, prim)
! var 2 --> Sx or vx
! var 3, --> Sy or vy

!assume same extent for both coordinates
	
!allocate(Unew(Nv,N,N),Uold(Nv,N,N),x(N,N),y(N,N),Pnew(Nv,N,N),Pold(Nv,N,N), &
!     &           vx_y(N,N),vy_x(N,N))!,fy(N,N),fx(N,N))
	
 	allocate(Unew(Nv,N,N),STAT=ierr)
    allocate(Uold(Nv,N,N),STAT=ierr)
    allocate(x(N,N),STAT=ierr)
    allocate(y(N,N),STAT=ierr)
    allocate(Pnew(Nv,N,N),STAT=ierr)
    allocate(Pold(Nv,N,N),STAT=ierr)
    allocate(vx_y(N,N),STAT=ierr)
    allocate(vy_x(N,N),STAT=ierr)!,fy(N,N),fx(N,N))

		
!define coords and related	
	
	dx= 10./(N-1)
	dy = dx
	dt = cfl * dx
	
	do i = 1, N
	x(i,:) =  (i-1)*dx
	y(:,i) =  (i-1)*dy
	end do
	
!define initial data
	call initial(Uold,Pold,dx,dy,x,y,Nv,N)

if (resume==1) time = 1.0d0*int(resume_time/dt)*dt
if (resume==0) time = 0.0d0
	Unew = Uold
	Pnew = Pold

	if(ltrace) print*, 'about to integrate'
	
	   do ITE=1, Nt
	    !integrate

        if(ltrace) print*, 'call evolve'
	    call evolve(unew,uold,Pnew,Pold,dx,dt,time,x,y,Nv,N,ITE,freq,arg)	
        if(ltrace) print*, 'back from evolve'

          if (mod(ITE-1,freq).eq.0) then
          print*, 'call out', time
!calculate vorticity
	  call derivs(pold(ivy,:,:),vy_x,dx,N,1)
	  call derivs(pold(ivx,:,:),vx_y,dx,N,2)
	  	  	
          ret = gft_out_brief(trim(arg)//'vort',time, (/n,n/), 2, vy_x-vx_y)  
          ret = gft_out_brief(trim(arg)//'vx',time, (/n,n/), 2, pold(ivx,:,:))
          ret = gft_out_brief(trim(arg)//'vy',time, (/n,n/), 2, pold(ivy,:,:))
          ret = gft_out_brief(trim(arg)//'eps',time, (/n,n/), 2, pold(irho,:,:))
!          ret = gft_out_brief('Sx',time, (/n,n/), 2, uold(ivx,:,:))
!          ret = gft_out_brief('Sy',time, (/n,n/), 2, uold(ivy,:,:))
!          ret = gft_out_brief('D',time, (/n,n/), 2, uold(irho,:,:))
 
        if(ltrace) print*, 'end call out'
	  end if

	   !update the dime and save the new values
	   !in the old array. 

            time = time + dt
            Uold = Unew
	    Pold = Pnew

	   end do

          ret = gft_out_brief(trim(arg)//'vx',time, (/n,n/), 2, pold(ivx,:,:))
          ret = gft_out_brief(trim(arg)//'vy',time, (/n,n/), 2, pold(ivy,:,:))
          ret = gft_out_brief(trim(arg)//'eps',time, (/n,n/), 2, pold(irho,:,:))
!          ret = gft_out_brief('Sx',time, (/n,n/), 2, uold(ivx,:,:))
!          ret = gft_out_brief('Sy',time, (/n,n/), 2, uold(ivy,:,:))
!          ret = gft_out_brief('D',time, (/n,n/), 2, uold(irho,:,:))

	end
