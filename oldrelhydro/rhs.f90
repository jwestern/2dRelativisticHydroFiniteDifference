 	MODULE M_RHS
	implicit none
	
	contains
	
	subroutine rhs(u,dxu,dyU,dxFxU,dyFyU,source,dx,Nx,x,pold,time,dt,printforcing,ITE,freq,arg)
	use m_pars, only : c, s,k_r, k_e, gamma,irho,ivx,ivy,solenoidalforce,rhoforce,fric,relativisticfriction, &
			& forcestr, fricstr, shiftrectforce, kforcel, kforceg, kfric, filterfric,lf,F0, &
			& gaussforce,shift,width, &
			& shiftgaussforce, shiftcosforce, shift2, width2, bofforce
	use fourier_tools
	use fourier_repackage
	use M_derivs

	implicit none
	real*8, dimension(:,:,:), intent(in) :: u, dxu,dyu, dxFxU,dyFyU,pold
	real*8, dimension(:,:,:), intent(inout) :: source
	real*8, dimension(:,:), intent(in) :: x
	real*8, intent(in) :: time, dt
	integer, intent(in) :: printforcing,ITE,freq
	character(len=*), intent(in) :: arg
	real*8 :: dx, v, ppi, e, wavenumber, alpha, beta, propfact
	integer :: i, j, Nx, ret, gft_out_brief
	logical :: ltrace

!My changes
	real*8, dimension(:,:), allocatable :: coordsx,coordsy,gmma,forcingx,forcingy,streamfunc,fricx,fricy
	real*8 :: extremevals(4)
	real(kind=8), dimension(:,:), allocatable :: storage,storage_x,storage_y

!Random seed getting/saving stuff
	integer :: seed_size
	integer, allocatable :: seed_master(:)

	seed_size=2

	allocate(coordsx(Nx,Nx),coordsy(Nx,Nx),gmma(Nx,Nx),forcingx(Nx,Nx), &
	&		 forcingy(Nx,Nx),streamfunc(Nx,Nx),fricx(Nx,Nx),fricy(Nx,Nx))

	if (mod(1.0d0*int(time/dt)*dt,1*dt).eq.0) then		!Skips every other time step.
!		write(*,*) 'first'
		call init_random_seed()			!create random seed from system clock
		call random_seed(size=seed_size)	!find out the size of the seed, save in seed_size
		allocate(seed_master(seed_size))	!allocate memory for the seed, to be saved in seed_master
		call random_seed(get=seed_master)	!get the seed, save in seed_master

		!Now save the seed to file, for later use.
		OPEN(1,FILE='seed.dat',STATUS='REPLACE',ACTION='WRITE')
		WRITE(1,*) seed_master
		CLOSE(1)
		deallocate(seed_master)
	else
!		write(*,*) 'else'
		call init_random_seed()			!generate random seed just to get its size
		call random_seed(size=seed_size)	!find out size and save (hope it's same size always)
		allocate(seed_master(seed_size))	!allocate
		
		!Read in old seed (i.e. overwrite seed_master)
		OPEN(1,FILE='seed.dat',ACTION='READ')
		READ(1,*) seed_master
		CLOSE(1)

		call random_seed(put=seed_master)	!re-initialize with old seed
		deallocate(seed_master)
	endif

!!!!!!!!!!!
	ppi=3.141592653589793
	e=2.718281828459045

!for debug
	ltrace=.true.
	

!!!!!!!! Do forcing function
!if (mod(int(time),1).eq.0) then
if (( bofforce==0.and.((gaussforce==0).and.(shiftgaussforce==0)) ).and.(shiftcosforce==0)) then
	if(solenoidalforce==0)then
		call forcingfield(forcingx,Nx,kforcel,kforceg,F0)
		call forcingfield(forcingy,Nx,kforcel,kforceg,F0)
	
!		extremevals=0.	
!		extremevals(1)=abs(maxval(forcingx))
!		extremevals(2)=abs(minval(forcingx))
!		forcingx=forcestr*.6*forcingx/maxval(extremevals(1:2))

!		extremevals=0.	
!		extremevals(1)=abs(maxval(forcingy))
!		extremevals(2)=abs(minval(forcingy))
!		forcingy=forcestr*.6*forcingy/maxval(extremevals(1:2))

		if(rhoforce==1)then
		forcingx=pold(irho,:,:)*forcingx
		forcingy=pold(irho,:,:)*forcingy
		endif

	elseif(solenoidalforce==1)then
		call forcingfield(streamfunc,Nx,kforcel,kforceg,F0)
		call derivs(streamfunc,forcingx,dx,Nx,2)
		call derivs(-streamfunc,forcingy,dx,Nx,1)

		!Collect extreme values for normalization
!		extremevals=0.
!		extremevals(1)=abs(minval(forcingx))
!		extremevals(2)=abs(maxval(forcingx))
!		extremevals(3)=abs(minval(forcingy))
!		extremevals(4)=abs(maxval(forcingy))
		
!		forcingx=forcestr*.6*forcingx/maxval(extremevals)
!		forcingy=-forcestr*.6*forcingy/maxval(extremevals)

		if(rhoforce==1)then
		forcingx=pold(irho,:,:)*forcingx
		forcingy=pold(irho,:,:)*forcingy
		endif
	endif

elseif (( bofforce==0.and.((gaussforce==1).and.(shiftgaussforce==0)) ).and.(shiftcosforce==0)) then

	if(solenoidalforce==0)then

		call forcingfield2(forcingx,Nx,lf,F0)
		call forcingfield2(forcingy,Nx,lf,F0)
	elseif(solenoidalforce==1)then

		call forcingfield2(streamfunc,Nx,lf,F0)
		call derivs(streamfunc,forcingx,dx,Nx,2)
                call derivs(-streamfunc,forcingy,dx,Nx,1)

        endif

elseif (( bofforce==0.and.((gaussforce==0).and.(shiftgaussforce==1)) ).and.(shiftcosforce==0)) then
	
	if(solenoidalforce==0)then

		call forcingfield3(forcingx,Nx,width,shift,F0)
		call forcingfield3(forcingy,Nx,width,shift,F0)
	elseif(solenoidalforce==1)then

                call forcingfield3(streamfunc,Nx,width,shift,F0)
                call derivs(streamfunc,forcingx,dx,Nx,2)
                call derivs(-streamfunc,forcingy,dx,Nx,1)

	endif

elseif (( bofforce==0.and.((gaussforce==0).and.(shiftgaussforce==0)) ).and.(shiftcosforce==1)) then

	if(solenoidalforce==0)then

		call forcingfield4(forcingx,Nx,width2,shift2,F0)
		call forcingfield4(forcingy,Nx,width2,shift2,F0)
	elseif(solenoidalforce==1)then

		call forcingfield4(streamfunc,Nx,width2,shift2,F0)
                call derivs(streamfunc,forcingx,dx,Nx,2)
                call derivs(-streamfunc,forcingy,dx,Nx,1)

        endif

elseif (( bofforce==1.and.((gaussforce==0).and.(shiftgaussforce==0)) ).and.(shiftcosforce==0)) then

	if(solenoidalforce==0)then

		call forcingfield5(forcingx,Nx,lf,F0)
		call forcingfield5(forcingy,Nx,lf,F0)
	elseif(solenoidalforce==1)then

		call forcingfield5(streamfunc,Nx,lf,F0)
		call derivs(streamfunc,forcingx,dx,Nx,2)
                call derivs(-streamfunc,forcingy,dx,Nx,1)

        endif

endif


!!!!!!!!!!!!! Make friction
if (fric==1) then
	if (relativisticfriction==1) then
		gmma=1./sqrt(1-pold(ivx,:,:)**2-pold(ivy,:,:)**2)

		fricx=(3./2.)*pold(irho,:,:)*gmma(:,:)**2*pold(ivx,:,:)
		fricy=(3./2.)*pold(irho,:,:)*gmma(:,:)**2*pold(ivy,:,:)
	elseif(relativisticfriction==0)then
		fricx=pold(ivx,:,:)
		fricy=pold(ivy,:,:)
	endif

	if (filterfric==1) then	!Apply Fourier filter if requested...
		call fftfilter(fricx,Nx,0,kfric)
		call fftfilter(fricy,Nx,0,kfric)
	endif
endif
!Dump fricx/y and forcingx/y to file, to check out...
if(1==0)then
        if (mod(ITE-1,freq).eq.0) then
	if(printforcing.eq.1)then
	ret = gft_out_brief('fricx',time, (/nx,nx/), 2, fricx)
	ret = gft_out_brief('fricy',time, (/nx,nx/), 2, fricy)
	endif
	endif
endif
        if (mod(ITE-1,freq).eq.0) then
	if(printforcing.eq.1)then
	ret = gft_out_brief(trim(arg)//'fx',time, (/nx,nx/), 2, forcingx)
	ret = gft_out_brief(trim(arg)//'fy',time, (/nx,nx/), 2, forcingy)
	endif
	endif

!Make edges match up
	forcingx(:,Nx)=forcingx(:,1)
	forcingx(Nx,:)=forcingx(1,:)
	forcingy(:,Nx)=forcingy(:,1)
	forcingy(Nx,:)=forcingy(1,:)
	fricx(:,Nx)=fricx(:,1)
	fricx(Nx,:)=fricx(1,:)
	fricy(:,Nx)=fricy(:,1)
	fricy(Nx,:)=fricy(1,:)



!!!!!!!!!!!
	source = 0.0
		
	source(irho,:,:) = -dxFxU(irho,:,:) - dyFyU(irho,:,:)

	source(ivx,:,:) = -dxFxU(ivx,:,:) - dyFyU(ivx,:,:) + forcingx - forcestr**(2./3.)*fricstr*fricx
	
	source(ivy,:,:) = -dxFxU(ivy,:,:) - dyFyU(ivy,:,:) + forcingy - forcestr**(2./3.)*fricstr*fricy

	deallocate(coordsx,coordsy,gmma,forcingx,forcingy,streamfunc,fricx,fricy)
			
	end subroutine rhs
	
	subroutine init_random_seed()
	INTEGER :: i, n, clock
	INTEGER, DIMENSION(:), ALLOCATABLE :: seed

	CALL RANDOM_SEED(size=n)
	ALLOCATE(seed(n))

	CALL SYSTEM_CLOCK(COUNT=clock)

	seed=clock+37*(/ (i-1, i=1,n) /)
	CALL RANDOM_SEED(PUT=seed)

	DEALLOCATE(seed)
	end subroutine

	INTEGER FUNCTION wrappedii(i,Nx)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: i	!coords of 'field point'.
		INTEGER :: N, newi
		INTEGER :: Nx
		N=Nx-1
		if (i<1) then
			newi=N+i
			wrappedii=newi
		elseif (i>N) then
			newi=i-N
			wrappedii=newi
		else
			wrappedii=i
		endif

	END FUNCTION wrappedii

	INTEGER FUNCTION wrappedjj(j,Nx)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: j	!coords of 'field point'.
		INTEGER :: N, newj
		INTEGER :: Nx
		N=Nx-1
		if (j<1) then
			newj=N+j
			wrappedjj=newj
		elseif (j>N) then
			newj=j-N
			wrappedjj=newj
		else
			wrappedjj=j
		endif
	END FUNCTION wrappedjj

	REAL FUNCTION wrappedi(x,i,N)
		IMPLICIT NONE
		REAL, INTENT(IN) :: x		!coords of peak.
		INTEGER, INTENT(IN) :: i	!coords of 'field point'.
		INTEGER :: N
		REAL :: gapx, newi
		N=N-1
if (1==0) then
		!Figure out how much displaced from middle.
		gapx=N-x			
		if ((N-x)<N/2) then		!In x...
			gapx=mod(gapx,float(N)/2)+1
		endif
endif
		!Find the wrapped coords of field point.
		if (abs(i-x)>200) then		!x...
			if ((N-x)<N/2) then
!				newi=-(N-i)
				newi=N+i
			else
!				newi=N+i
				newi=-(N-i)
			endif
			wrappedi=newi
		else
			wrappedi=i
		endif		
	END FUNCTION wrappedi

	REAL FUNCTION wrappedj(y,j,N)
		IMPLICIT NONE
		REAL, INTENT(IN) :: y		!coords of peak.
		INTEGER, INTENT(IN) :: j	!coords of 'field point'.
		INTEGER :: N
		REAL :: gapy, newj
		N=N-1
if (1==0) then
		!Figure out how much displaced from middle.
		gapy=N-y			!In y...
		if ((N-y)<N/2) then
			gapy=mod(gapy,float(N)/2)+1
		endif
endif
		!Find the wrapped coords of field point.
		if (abs(j-y)>200) then		!y...
			if ((N-y)<N/2) then
!				newj=-(N-j)
				newj=N+j
			else
				!newj=N+j
				newj=-(N-j)
			endif
			wrappedj=newj
		else
			wrappedj=j
		endif
	END FUNCTION wrappedj

	subroutine bumpp(bump,Nx)
		REAL, DIMENSION(0:2) :: coords
		INTEGER :: Nx
		REAL, DIMENSION(:,:),allocatable :: bump
		REAL, PARAMETER :: e=2.718281828459045
		INTEGER :: i,j
		allocate(bump(Nx,Nx))
		call init_random_seed()
		call random_number(coords)
		coords=coords*400.
		
		do i=1,Nx
			do j=1,Nx
				bump(i,j)=3.*(10.**(-1))*e**(-(coords(1)-wrappedi(coords(1),i,Nx))**2/(2.*(3.)**2)-(coords(2)-wrappedj(coords(2),j,Nx))**2/(2.*(3.)**2))
			enddo
		enddo
		deallocate(bump)
	end subroutine

	REAL FUNCTION G(kx,ky)
		IMPLICIT NONE
		INTEGER, INTENT(IN) :: kx, ky
		if ((sqrt(float(kx**2+ky**2))<191).and.(sqrt(float(kx**2+ky**2))>143)) then
			G=1.
		else
			G=0.
		endif
	END FUNCTION G

	end module M_RHS
