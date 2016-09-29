	MODULE M_PARS
	implicit none
	
!define parameters !!!!!!!!!!!!!!!!!!!!!!!
	real*8 ::  cfl, weos, disip, sL, sR, c,s,a,k_r,k_e, gamma, cl,valuefloor,fricstr,resume_time,F0
	integer :: N, Nt, freq, bc, derorder,dissorder,STMETH,solenoidalforce,rhoforce, &
			&fric,relativisticfriction, forcestr, shiftrectforce, kforcel, kforceg, resume, &
			&resume_integer, kfric, filterfric, lf,gaussforce,shift,width,shiftgaussforce, &
			&shiftcosforce, shift2, width2, bofforce
	integer :: irho, ivx,ivy	
	
	CONTAINS
	
	
	subroutine readpars
	implicit none
        namelist /pars_input/ N, Nt, freq, cfl, weos, disip, &
	&		sL, sR, derorder,dissorder, c,s,a, STMETH, &
	&		k_r, k_e, gamma, cl,valuefloor, &
	&		rhoforce, fric, relativisticfriction, forcestr, fricstr, &
	&		kfric, filterfric, resume, resume_integer, resume_time, &
	&		F0, solenoidalforce, shiftrectforce, kforcel, kforceg, &
	&		gaussforce, lf, shiftgaussforce, shift, width, &
	&		shiftcosforce, shift2, width2, bofforce

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
