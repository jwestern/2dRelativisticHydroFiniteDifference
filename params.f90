	MODULE M_PARS
	implicit none
	
!define parameters !!!!!!!!!!!!!!!!!!!!!!!
	real*8 ::  cfl, weos, disip, sL, sR, c,s,a,k_r,k_e, gamma, cl,valuefloor,fricstr,&
		&	resume_time,F0,Lx
	integer :: N, Nt, freq, derorder,dissorder,STMETH,solenoidalforce, &
		&	fric,relativisticfriction, shiftrectforce, kforcel, kforceg, resume, &
		&	resume_integer, kfric, filterfric, lf,gaussforce, printforcing,&
		&	seed
	integer(kind=8) :: nskip
	integer :: irho, ivx,ivy	
	
	CONTAINS
	
	
	subroutine readpars
	implicit none
        namelist /pars_input/ Lx, N, Nt, freq, cfl, weos, disip, &
	&		sL, sR, derorder,dissorder, STMETH,c,s,a, &
	&		k_r, k_e, gamma, cl,valuefloor, &
	&		fric, relativisticfriction, fricstr, &
	&		kfric, filterfric, resume, resume_integer, resume_time, &
	&		printforcing, seed, nskip, &
	&		F0, solenoidalforce, shiftrectforce, kforcel, kforceg, &
	&		gaussforce, lf

!read params !!!!!!!!!!!!!!!!!!!!!!!


   open (unit = 10, file = "pars.in", status = "old" )
   read (unit = 10, nml = pars_input)
   close(unit = 10)

!	open (unit=20, file = "pars2.in", status = "replace")
!	write (20,*) N
!	close (unit=20)

	irho = 1
	ivx = 2
	ivy = 3

	
	end subroutine readpars
	
	
	end module M_pars	
