	MODULE M_DERIVS
	
	implicit none
	
	
	contains
	
!this is just a wrapper, it calls the ders later	
	subroutine derivs(u,du,dx, Nx,dir)
	use m_pars, only : derorder
	real*8, dimension(:,:), intent(in) :: u
	real*8, dimension(:,:), intent(inout) :: du
	real*8 :: dx
	integer :: i, j, nx,dir
	
	
	if(dir.eq.1) then
	    do j = 1, Nx
		if (derorder.eq.1) then
		  call derivs12(u(:,j),du(:,j),dx,NX)
		else if(derorder.eq.2) then
		  call derivs24(u(:,j),du(:,j),dx,NX)
		else if(derorder.eq.3) then
		  call derivs36(u(:,j),du(:,j),dx,NX)
		else if(derorder.eq.4) then
		  call derivs48(u(:,j),du(:,j),dx,NX)
		else if(derorder.eq.41) then
		  call derivs48_b(u(:,j),du(:,j),dx,NX)
		end if
	    end do	
	else 
	
	    do j = 1, Nx
		if (derorder.eq.1) then
		  call derivs12(u(j,:),du(j,:),dx,NX)
		else if(derorder.eq.2) then
		  call derivs24(u(j,:),du(j,:),dx,NX)
		else if(derorder.eq.3) then
		  call derivs36(u(j,:),du(j,:),dx,NX)
		else if(derorder.eq.4) then
		  call derivs48(u(j,:),du(j,:),dx,NX)
		else if(derorder.eq.41) then
		  call derivs48_b(u(j,:),du(j,:),dx,NX)
		end if
	    end do
	    	
	end if
		
	end subroutine derivs
	
	subroutine derivs12(u,du,dx, Nx)
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, nx
	
	do i = 2, nx-1
	du(i) = (u(i+1)-u(i-1))/dx*0.5
	end do
	du(1) = (u(2)-u(1))/dx
	du(nx) = (u(nx)-u(nx-1))/dx
	
!use periodic bcs
	du(1) = 0.5*(u(2)-u(nx-1))/dx
	du(nx) = 0.5*(u(2)-u(nx-1))/dx

	end subroutine derivs12
	
		
	subroutine dissip(u,du,dx, Nx)
	use m_pars, only : dissorder
	
	real*8, dimension(:,:), intent(in) :: u
	real*8, dimension(:,:), intent(inout) :: du
	real*8 :: dx
	integer :: i, j, Nx
	
	real*8, allocatable, dimension(:,:) :: tempu
	
	allocate(tempu(Nx,Nx))
	
!do each dimension and then add
	do j=1, Nx
		if (dissorder.eq.4) then
		    call dissip4(u(:,j),du(:,j),dx,NX)
		else if(dissorder.eq.6) then
		    call dissip6(u(:,j),du(:,j),dx,NX)
		else if(dissorder.eq.8) then
		   call dissip8(u(:,j),du(:,j),dx,NX)	
		end if
	end do
	
	do j =1, Nx
		if (dissorder.eq.4) then
		    call dissip4(u(j,:),tempu(j,:),dx,NX)
		else if(dissorder.eq.6) then
		    call dissip6(u(j,:),tempu(j,:),dx,NX)
		else if(dissorder.eq.8) then
		   call dissip8(u(j,:),tempu(j,:),dx,NX)	
		end if	
	end do
	
	du(1:nx,1:nx) = (du(1:nx,1:nx)+tempu(1:nx,1:nx))/dx
	
	deallocate(tempu)
	
	end subroutine dissip
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!		DISIPATION FIRST 4th order interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine dissip4(u,du,dx, Nx)
	use m_pars, only : derorder		
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, Nx
	real*8 :: smr3,smr2,smr1,spr1,spr2,spr3
	
				
	du = 0.0
		
        do i = 3, nx-2
        du(i) = u(i+2)+u(i-2)-4.*(u(i+1)+u(i-1))+6*u(i)
        end do
        du(2) = u(4)-4.*u(3)+5.*u(2)-2.*u(1)
        du(1) = 2.*(u(3)-2.*u(2)+u(1))
        du(nx-1) = u(nx-3)-4.*u(nx-2)+5.*u(nx-1)-2.*u(nx)
        du(nx) = 2.*(u(nx)-2.*u(nx-1)+u(nx-2))

!force periodicity
	du(2) = u(4)+u(nx-1)-4.*(u(3)+u(1))+6*u(2)
        du(1) = u(3)+u(nx-2)-4.*(u(2)+u(nx-1))+6*u(1)
	du(nx) = du(1)
        du(nx-1) = u(2)+u(nx-3)-4.*(u(1)+u(nx-2))+6*u(nx-1)


	
	end subroutine dissip4	
	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   DISIPATION 6th order interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine dissip6(u,du,dx, Nx)	
	use m_pars, only : derorder			
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, Nx, lef, rig			
	real*8 :: smr3,smr2,smr1,spr1,spr2,spr3

!fix scalar product according to derivative operators
	if(derorder.eq.1) then
	lef=4; rig=nx-3
	smr3=0.5;smr2=1.;smr1=1.
	spr3=smr3;spr2=smr2;spr1=smr1
	else if (derorder.eq.2) then
	lef=5; rig=nx-4	
	smr3=59.0/48.0;smr2=43.0/48.0;smr1=49.0/48.0
	spr3=smr3;spr2=smr2;spr1=smr1
	else if (derorder.eq.3) then
	lef=7;rig=nx-6
	smr3=5359.0/4320.0;smr2=7877.0/8640.0;smr1=43801.0/43200.0
	spr3=smr3;spr2=smr2;spr1=smr1	
	else if (derorder.eq.41) then
	lef=9;rig=nx-8
	smr3=103097.0/80640.0;smr2=670091.0/725760.0;smr1=5127739.0/5080320.0
	spr3=smr3;spr2=smr2;spr1=smr1
	end if
			
	du = 0.0	

!now get the dissipation
        do i = lef, rig
        du(i) = u(i+3)+u(i-3)-6.*(u(i+2)+u(i-2))+ 15.*(u(i+1)+u(i-1))-20.*u(i)
        end do
	
        du(lef-3) = (u(lef)-3.*u(lef-1)+3.*u(lef-2)-u(lef-3))/smr3
        du(lef-2) = (u(lef+1)-6.*u(lef)+12.*u(lef-1)-10.*u(lef-2)+3.*u(lef-3))/smr2
        du(lef-1) = (u(lef+2)-6.*u(lef+1)+15.*u(lef)-19.*u(lef-1)+12.*u(lef-2)-3.*u(lef-3))/smr1

        du(rig+3) = (u(rig)-3.*u(rig+1)+3.*u(rig+2)-u(rig+3))/spr3
        du(rig+2) = (u(rig-1)-6.*u(rig)+12.*u(rig+1)-10.*u(rig+2)+3.*u(rig+3))/spr2
        du(rig+1) = (u(rig-2)-6.*u(rig-1)+15.*u(rig)-19.*u(rig+1)+12.*u(rig+2)-3.*u(rig+3))/spr1

	
	end subroutine dissip6	
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   DISIPATION 8th order interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine dissip8(u,du,dx, Nx)	
	use m_pars, only : derorder			
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, Nx, lef, rig			
	real*8 :: smr4,smr3,smr2,smr1,spr1,spr2,spr3,spr4

!fix scalar product according to derivative operators
	if(derorder.eq.1) then
	lef=5; rig=nx-4
	smr4=0.5;smr3=1.;smr2=1.;smr1=1.
	spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1
	else if (derorder.eq.2) then
	lef=5; rig=nx-4	
	smr4=17./48.;smr3=59.0/48.0;smr2=43.0/48.0;smr1=49.0/48.0
	spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1
	else if (derorder.eq.3) then
	lef=7;rig=nx-6
	smr4=2711.0/4320.0;smr3=5359.0/4320.0;smr2=7877.0/8640.0;smr1=43801.0/43200.0
	spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1	
	else if (derorder.eq.41) then
	lef=9;rig=nx-8
	smr4 = 299527.0/725760.0;smr3=103097.0/80640.0;smr2=670091.0/725760.0;smr1=5127739.0/5080320.0
	spr4=smr4;spr3=smr3;spr2=smr2;spr1=smr1
	end if
			
	du = 0.0	

!now get the dissipation
        do i = lef, rig
        du(i) = -(u(i+4)+u(i-4)-8.*(u(i+3)+u(i-3))+28.*(u(i+2)+u(i-2)) &
	&	-56.*(u(i+1)+u(i-1))+70.*u(i))
        end do
	
        du(lef-4) = (-u(lef)+4.*u(lef-1)-6.*u(lef-2)+4.*u(lef-3)-u(lef-4))/smr4
        du(lef-3) = (2.*u(lef)-9*u(lef-1)+15.*u(lef-2)-11.*u(lef-3)+3.*u(lef-4))/smr3
        du(lef-2) = (-u(lef+1)+3.*u(lef)-8.*u(lef-2)+9.*u(lef-3)-3.*u(lef-4))/smr2
        du(lef-1) = (-u(lef+2)+6.*u(lef+1)-14.*u(lef)+15.*u(lef-1)-6.*u(lef-2)-u(lef-3)+u(lef-4))/smr1

        du(rig+4) = (-u(rig)+4.*u(rig+1)-6.*u(rig+2)+4.*u(rig+3)-u(rig+4))/spr4
        du(rig+3) = (2.*u(rig)-9*u(rig+1)+15.*u(rig+2)-11.*u(rig+3)+3.*u(rig+4))/spr3
        du(rig+2) = (-u(rig-1)+3.*u(rig)-8.*u(rig+2)+9.*u(rig+3)-3.*u(rig+4))/spr2
        du(rig+1) = (-u(rig-2)+6.*u(rig-1)-14.*u(rig)+15.*u(rig+1)-6.*u(rig+2)-u(rig+3)+u(rig+4))/spr1
	
	end subroutine dissip8		
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!		DERIVATIVES FROM NOW ON	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! 2nd in the bound, 4th in the interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine derivs24(u,du,dx,Nx)
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, nx

	real*8 :: d00,d01,d02,d03,&
		 & d10, d12, &
		 & d20,d21,d23,d24, &
		 & d30,d32,d34,d35, &
		 & d1, d2
		 

   d1 = 0.6666666666666666666666
   d2 = -0.083333333333333333333

!  d00 = -1.4117647058823529411764705882353
!  d01 = 1.7352941176470588235294117647059
!  d02 = -0.23529411764705882352941176470588
!  d03 = -0.088235294117647058823529411764706
!
! d10 = -0.50000000000000000000000000000000
! d12 = 0.50000000000000000000000000000000
!
! d20 = 0.093023255813953488372093023255814
! d21 = -0.68604651162790697674418604651163
! d23 = 0.68604651162790697674418604651163
! d24 = -0.093023255813953488372093023255814
!
!! d30 = 0.030612244897959183673469387755102
! d32 = -0.60204081632653061224489795918367
! d34 = 0.65306122448979591836734693877551
! d35 = -0.081632653061224489795918367346939
!	
!	du(1) = (d00*u(1)+d01*u(2)+d02*u(3)+d03*u(4))/dx
!	du(2) = (d10*u(1)+d12*u(3))/dx
!	du(3) = (d20*u(1)+d21*u(2)+d23*u(4)+d24*u(5))/dx
!	du(4) = (d30*u(1)+d32*u(3)+d34*u(5)+d35*u(6))/dx
	
	
	do i = 3, nx-2
	du(i) = (-d2*u(i-2)-d1*u(i-1)+d1*u(i+1)+d2*u(i+2))/dx
	end do

        du(1) =  (-d2*u(nx-2)-d1*u(nx-1)+d1*u(2)+d2*u(3))/dx
        du(2) =  (-d2*u(nx-1)-d1*u(nx)+d1*u(3)+d2*u(4))/dx
        du(nx) = (-d2*u(nx-2)-d1*u(nx-1)+d1*u(2)+d2*u(3))/dx
        du(nx-1) = (-d2*u(nx-3)-d1*u(nx-2)+d1*u(1)+d2*u(2))/dx

!	du(nx) = -(d00*u(nx)+d01*u(nx-1)+d02*u(nx-2)+d03*u(nx-3))/dx
!	du(nx-1) = -(d10*u(nx)+d12*u(nx-2))/dx
!	du(nx-2) = -(d20*u(nx)+d21*u(nx-1)+d23*u(nx-3)+d24*u(nx-4))/dx
!	du(nx-3) = -(d30*u(nx)+d32*u(nx-2)+d34*u(nx-4)+d35*u(nx-5))/dx
			
	end subroutine derivs24

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! 3rth in the bound, 6th in the interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	subroutine derivs36(u,du,dx,Nx)
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, nx

	real*8 :: d00,d01,d02,d03,d04, &
		 & d10, d12, d13,d14,d15,&
		 & d20,d21,d23,d24,d25, &
		 & d30,d31,d32,d34,d35,d36, &
		 & d40,d41,d42,d43,d44,d45,d46,d47,&
		 & d50,d51,d52,d53,d54,d55,d56,d57,d58,&
		 & d1, d2, d3
		 
 
d00= -1.5825335189391164187852589933328
d01= 1.9968007424231323418077026399980
d02= 0.0047988863653014872884460400029306
d03= -0.66986592424353432485896402666862
d04= 0.25079981439421691454807434000049

d10= -0.45374732928216654180193679069897
d12= 0.20413995948833208468603457365632
d13= 0.42505341435666916396126418602070
d14= -0.19379006076750187297094813951552
d15= 0.018344016204667166125586170537473

d20= -0.0024160826263371449649575802286979
d21= -0.45229312676749047092093938276159
d23= 0.23791958686831427517521209885651
d24= 0.34541374646501905815812123447682
d25= -0.12862412393950571744743637034305

d30= 0.17061018846799776077626422840082
d31= -0.47641039995023947253840890713442
d32= -0.12035827579772345586863220750140
d34= 0.42710082726876904895191889034024
d35= -0.014377682403433476394849785407725
d36= 0.013435342414629595073707781302482

d40= -0.086915492361728238331005882104016
d41= 0.29554398882823409927637425415767
d42= -0.23775972239854428504929964876645
d43= -0.58114341331302103169565401379544
d45= 0.75652321103635055647243028225636
d46= -0.16452964326520248825695061571664
d47= 0.018281071473911387584105623968516

d51= -0.025155437851495019139593464380570
d52= 0.079610054564964270222141047008059
d53= 0.017590922581676217437958037487729
d54= -0.68025083141176381056749084876297
d56= 0.73970913906075203762471176457159
d57= -0.14794182781215040752494235291432
d58= 0.016437980868016711947215816990480



d1= 0.750000000000000000000000000000000
d2= -0.15000000000000000000000000000000
d3= 0.016666666666666666666666666666667  
	
	du(1) = (d00*u(1)+d01*u(2)+d02*u(3)+d03*u(4)+d04*u(5))/dx
	du(2) = (d10*u(1)+d12*u(3)+d13*u(4)+d14*u(5)+d15*u(6))/dx
	du(3) = (d20*u(1)+d21*u(2)+d23*u(4)+d24*u(5)+d25*u(6))/dx
	du(4) = (d30*u(1)+d31*u(2)+d32*u(3)+d34*u(5)+d35*u(6)+d36*u(7))/dx
	du(5) = (d40*u(1)+d41*u(2)+d42*u(3)+d43*u(4)+d45*u(6)+d46*u(7)+d47*u(8))/dx
	du(6) = (d51*u(2)+d52*u(3)+d53*u(4)+d54*u(5)+d56*u(7)+d57*u(8)+d58*u(9))/dx
	
	
	do i = 7, nx-6
	du(i) = (-d3*u(i-3)-d2*u(i-2)-d1*u(i-1)+d1*u(i+1)+d2*u(i+2)+d3*u(i+3))/dx
	end do
	

	du(nx) = -(d00*u(nx)+d01*u(nx-1)+d02*u(nx-2)+d03*u(nx-3)+d04*u(nx-4))/dx
	du(nx-1) = -(d10*u(nx)+d12*u(nx-2)+d13*u(nx-3)+d14*u(nx-4)+d15*u(nx-5))/dx
	du(nx-2) = -(d20*u(nx)+d21*u(nx-1)+d23*u(nx-3)+d24*u(nx-4)+d25*u(nx-5))/dx
	du(nx-3) = -(d30*u(nx)+d31*u(nx-1)+d32*u(nx-2)+d34*u(nx-4)+d35*u(nx-5)+d36*u(nx-6))/dx
	du(nx-4) = -(d40*u(nx)+d41*u(nx-1)+d42*u(nx-2)+d43*u(nx-3)+d45*u(nx-5)+d46*u(nx-6)+d47*u(nx-7))/dx
	du(nx-5) = -(d51*u(nx-1)+d52*u(nx-2)+d53*u(nx-3)+d54*u(nx-4)+d56*u(nx-6)+d57*u(nx-7)+d58*u(nx-8))/dx
			
	end subroutine derivs36
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! 4th in the bound, 8th in the interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	subroutine derivs48(u,du,dx,Nx)
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, nx		

   real*8 :: d00,d01,d02,d03,d04,d05,d06,d07,d08,d09,d010,d011, &
        &    d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d110,d111, &
        &    d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d210,d211, &
        &    d30,d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,d311, &
        &    d40,d41,d42,d43,d44,d45,d46,d47,d48,d49,d410,d411, &
        &    d50,d51,d52,d53,d54,d55,d56,d57,d58,d59,d510,d511, &
        &    d60,d61,d62,d63,d64,d65,d66,d67,d68,d69,d610,d611, &
        &    d70,d71,d72,d73,d74,d75,d76,d77,d78,d79,d710,d711, dc1,dc2,dc3,dc4


   d00=-1.6955436044318985087498556542484
   d01=2.0610513554928258770826116045752
   d02=0.87789728901434824583477679084963
   d03=-2.5445639556810149125014434575163
   d04=1.6889486445071741229173883954248
   d05=-0.38778972890143482458347767908496

   d10=-0.39835918734972926809517745906661
   d12=-0.44127885642072764523900478066757
   d13=1.8989658393441626095262349706691
   d14=-1.5738365692621829357170143420027
   d15=0.60604617027316423238240764906811
   d16=-0.091537396584686992857446038000302

   d20=-1.0055577090642838920860210778808
   d21=2.6151125596961539353524325346492
   d23=-10.448835244859039473942416960576
   d24=16.854276269429682078398421541815
   d25=-10.479793844227156688020808246231
   d26=2.5403284648744886691035720465902
   d27=-0.075530495849844628805179838367064

   d30=0.41730852995727528198434306408336
   d31=-1.6112948490126556929673253900322
   d32=1.4960581706703734383449080548465
   d34=-1.6738166082885887268879133240715
   d35=1.9668117242490862700533359284822
   d36=-0.64585509260926636725392126737262
   d37=0.050788125033775796726572934064281

   d40=-1.2067978767157806109219504339143
   d41=5.8182409265274916785465083281307
   d42=-10.513925845304986417473772536922
   d43=7.2925946575767793888363987219850
   d45=-3.0385143353798934097204370001146
   d46=2.0057652009112144592418490931814
   d47=-0.34870908370292399111170040001164
   d48=-0.0086536439119010973968957723343805
  
   d50=0.089446187550400088107295882497877
   d51=-0.72324463225750363106445242685870
   d52=2.1103523865873885758072494834961
   d53=-2.7662054325681783318765966175687
   d54=0.98086385026111691291085453864196
   d56=0.30802293704705996043280348352271
   d57=-0.026238992169135626268229544755535
   d58=0.029797181295285022842565739061272
   d59=-0.0027934857464329708914905380369943

   d61=0.15126303740835199995224529205735
   d62=-0.70834831886017471258878769201994
   d63=1.2577996869081960509841200672744
   d64=-0.89656603854302375846464634405377
   d65=-0.42651843804299217071512177699248
   d67=0.80159349003841766765011518336067
   d68=-0.21661535522787203529072916962025
   d69=0.041260067662451816245853175165761
   d610=-0.0038681313433548577730487351717901

   d72=0.019265719907611002302054245220619
   d73=-0.090478311526126167237971095382715
   d74=0.14258418768974005892265577479665
   d75=0.033235928479719164073418453370319
   d76=-0.73326354624003548802568409455577
   d78=0.79260196355547737511601116983528
   d79=-0.19815049088886934377900279245882
   d710=0.037742950645498922624571960468347
   d711=-0.0035384016230155239960536212939075
  
   dc1= 0.800000000000000000000000000000000
   dc2= -0.2000000000000000000000000000000
   dc3= 0.038095238095238095238095238095238
   dc4= -0.003571428571428571428571428571428
  
!now get busy

  du(1) = (d00*u(1)+d01*u(2)+d02*u(3)+d03*u(4)+d04*u(5)+d05*u(6))/dx
  du(2) = (d10*u(1)+d12*u(3)+d13*u(4)+d14*u(5)+d15*u(6)+d16*u(7))/dx
  du(3) = (d20*u(1)+d21*u(2)+d23*u(4)+d24*u(5)+d25*u(6)+d26*u(7)+d27*u(8))/dx
  du(4) = (d30*u(1)+d31*u(2)+d32*u(3)+d34*u(5)+d35*u(6)+d36*u(7)+d37*u(8))/dx		
  du(5) = (d40*u(1)+d41*u(2)+d42*u(3)+d43*u(4)+d45*u(6)+d46*u(7)+d47*u(8)+d48*u(9))/dx		
  du(6) = (d50*u(1)+d51*u(2)+d52*u(3)+d53*u(4)+d54*u(5)+d56*u(7)+d57*u(8)+d58*u(9)+d59*u(10))/dx		
  du(7) = (d61*u(2)+d62*u(3)+d63*u(4)+d64*u(5)+d65*u(6)+d67*u(8)+d68*u(9)+d69*u(10)+d610*u(11))/dx				
  du(8) = (d72*u(3)+d73*u(4)+d74*u(5)+d75*u(6)+d76*u(7)+d78*u(9)+d79*u(10)+d710*u(11)+d711*u(12))/dx		
	
  do i=9,nx-8
  du(i) = (-dc4*u(i-4)-dc3*u(i-3)-dc2*u(i-2)-dc1*u(i-1) &
  &	   +dc4*u(i+4)+dc3*u(i+3)+dc2*u(i+2)+dc1*u(i+1))/dx
  end do		
	

  du(nx)   = -(d00*u(nx)+d01*u(nx-1)+d02*u(nx-2)+d03*u(nx-3)+d04*u(nx-4)+d05*u(nx-5))/dx
  du(nx-1) = -(d10*u(nx)+d12*u(nx-2)+d13*u(nx-3)+d14*u(nx-4)+d15*u(nx-5)+d16*u(nx-6))/dx
  du(nx-2) = -(d20*u(nx)+d21*u(nx-1)+d23*u(nx-3)+d24*u(nx-4)+d25*u(nx-5)+d26*u(nx-6)+d27*u(nx-7))/dx
  du(nx-3) = -(d30*u(nx)+d31*u(nx-1)+d32*u(nx-2)+d34*u(nx-4)+d35*u(nx-5)+d36*u(nx-6)+d37*u(nx-7))/dx		
  du(nx-4) = -(d40*u(nx)+d41*u(nx-1)+d42*u(nx-2)+d43*u(nx-3)+d45*u(nx-5)+d46*u(nx-6)+d47*u(nx-7)+d48*u(nx-8))/dx		
  du(nx-5) = -(d50*u(nx)+d51*u(nx-1)+d52*u(nx-2)+d53*u(nx-3)+d54*u(nx-4)+d56*u(nx-6)+d57*u(nx-7)+d58*u(nx-8)+d59*u(nx-9))/dx		
  du(nx-6) = -(d61*u(nx-1)+d62*u(nx-2)+d63*u(nx-3)+d64*u(nx-4)+d65*u(nx-5)+d67*u(nx-7)+d68*u(nx-8)+d69*u(nx-9)+d610*u(nx-10))/dx
  du(nx-7) = -(d72*u(nx-2)+d73*u(nx-3)+d74*u(nx-4)+d75*u(nx-5)+d76*u(nx-6)+d78*u(nx-8)+d79*u(nx-9)+d710*u(nx-10)+d711*u(nx-11))/dx
  		
	end subroutine derivs48
		
		
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!! 4th in the bound, 8th in the interior
!!!!!!!!!!!!!!!!!!!!!!!!!!!		
	subroutine derivs48_b(u,du,dx,Nx)
	real*8, dimension(:), intent(in) :: u
	real*8, dimension(:), intent(inout) :: du
	real*8 :: dx
	integer :: i, nx		

   real*8 :: d00,d01,d02,d03,d04,d05,d06,d07,d08,d09,d010,d011, &
        &    d10,d11,d12,d13,d14,d15,d16,d17,d18,d19,d110,d111, &
        &    d20,d21,d22,d23,d24,d25,d26,d27,d28,d29,d210,d211, &
        &    d30,d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,d311, &
        &    d40,d41,d42,d43,d44,d45,d46,d47,d48,d49,d410,d411, &
        &    d50,d51,d52,d53,d54,d55,d56,d57,d58,d59,d510,d511, &
        &    d60,d61,d62,d63,d64,d65,d66,d67,d68,d69,d610,d611, &
        &    d70,d71,d72,d73,d74,d75,d76,d77,d78,d79,d710,d711, dc1,dc2,dc3,dc4

   d00= -1.6955436044318985087498556542484
   d01= 2.252142137678813514633822362278
   d02= -0.05888015064022764242837280119
   d03= -0.72703849019795003445385686285
   d04= -0.03519446459907925766567721687
   d05= 0.3808994692748803682435341446956
   d06= -0.097708425809176140086689775338
   d07= -0.018676471275362299492904196473
   d08= 0.
   d09= 0.
   d010=0.
   d011= 0.

   d10= -0.4352931378302752275823881078793
   d11= 0.
   d12= 0.099854313212144418846805809049
   d13= 0.48598825799891087114955473054
   d14= -0.04056967339078804101602657017
   d15= -0.1516077655067655130871564977011
   d16= 0.028751917941456163466861493696
   d17= 0.012876087575317328222349142474
   d18= 0.
   d19= 0.
   d110= 0.
   d111= 0.

   d20= 0.06744227386055814358882611857
   d21= -0.59175794358011014241446301559
   d22= 0.
   d23= 0.09922290191545044394136441726
   d24= 1.24445434548753271358155515947
   d25= -1.368619045325369683541255238187
   d26= 0.68543519098309329993738259237
   d27= -0.13617772334115477509341003392
   d28= 0.
   d29= 0.
   d210= 0.
   d211= 0.

   d30= 0.119234324171531019965386581760
   d31= -0.41236675277145806241997120304
   d32= -0.01420667755300961916445912985
   d33= 0.
   d34= -0.117132848377661592077848335517
   d35= 0.6750458046328004849045861337825
   d36= -0.288099428939536862991210465807
   d37= 0.037525578837334631783516418677
   d38= 0.
   d39= 0.
   d310= 0.
   d311= 0.

   d40= 0.025147363295176347088194768037
   d41= 0.14998007970344354042651691946
   d42= -0.77630747812384192410033152270
   d43= 0.51033212364828546341398271274
   d44= 0.
   d45= -0.568060141489748837333529197702
   d46= 0.891763513806768671939424492616
   d47= -0.224201816928182164037362400114
   d48= -0.0086536439119010973968957723343805
   d49= 0.
   d410= 0.
   d411= 0.

   d50= -0.0878569049859194092294958469532
   d51= 0.1809259887937250033140311227933
   d52= 0.275603557814485387547649301144
   d53= -0.949412365700909499468138419805
   d54= 0.1833756882676831850910631088522
   d55= 0.
   d56= 0.42315722086966643064298670184389
   d57= -0.052796880607583149849171168899192
   d58= 0.029797181295285022842565739061272
   d59= -0.0027934857464329708914905380369943
   d510= 0.
   d511= 0.

   d60= 0.031207020141045519688121961545
   d61= -0.047511755865994320174424070761
   d62= -0.191127593117949651614482212113
   d63= 0.561072252375672358928364455971
   d64= -0.398613397284846386535560095569
   d65= -0.58594453589139385546142240382277
   d66= 0.
   d67= 0.81014142855224141198732709437972
   d68= -0.21661535522787203529072916962025
   d69= 0.041260067662451816245853175165761
   d610= -0.0038681313433548577730487351717901
   d611= 0.

   d70= 0.0054565862264050490869367571164
   d71= -0.019463641447689387726897436343
   d72= 0.034735133749982204632490070184
   d73= -0.066851276946818081029475174146
   d74= 0.091674221978406727279476067978
   d75= 0.066875790674993403525413442454852
   d76= -0.74108283592437134573347044379599
   d77= 0.
   d78= 0.79260196355547737511601116983528
   d79= -0.19815049088886934377900279245882
   d710= 0.037742950645498922624571960468347
   d711= -0.0035384016230155239960536212939075

   dc1= 0.800000000000000000000000000000000
   dc2= -0.2000000000000000000000000000000
   dc3= 0.038095238095238095238095238095238
   dc4= -0.0035714285714285714285714285714286 
!now get busy

  du(1) = (d00*u(1)+d01*u(2)+d02*u(3)+d03*u(4)+d04*u(5)+d05*u(6)+d06*u(7)+d07*u(8)+d08*u(9)+d09*u(10)+d010*u(11)+d011*u(12))/dx
  du(2) = (d10*u(1)+d11*u(2)+d12*u(3)+d13*u(4)+d14*u(5)+d15*u(6)+d16*u(7)+d17*u(8)+d18*u(9)+d19*u(10)+d110*u(11)+d111*u(12))/dx
  du(3) = (d20*u(1)+d21*u(2)+d22*u(3)+d23*u(4)+d24*u(5)+d25*u(6)+d26*u(7)+d27*u(8)+d28*u(9)+d29*u(10)+d210*u(11)+d211*u(12))/dx
  du(4) = (d30*u(1)+d31*u(2)+d32*u(3)+d33*u(4)+d34*u(5)+d35*u(6)+d36*u(7)+d37*u(8)+d38*u(9)+d39*u(10)+d310*u(11)+d311*u(12))/dx		
  du(5) = (d40*u(1)+d41*u(2)+d42*u(3)+d43*u(4)+d44*u(5)+d45*u(6)+d46*u(7)+d47*u(8)+d48*u(9)+d49*u(10)+d410*u(11)+d411*u(12))/dx		
  du(6) = (d50*u(1)+d51*u(2)+d52*u(3)+d53*u(4)+d54*u(5)+d55*u(6)+d56*u(7)+d57*u(8)+d58*u(9)+d59*u(10)+d510*u(11)+d511*u(12))/dx		
  du(7) = (d60*u(1)+d61*u(2)+d62*u(3)+d63*u(4)+d64*u(5)+d65*u(6)+d66*u(7)+d67*u(8)+d68*u(9)+d69*u(10)+d610*u(11)+d611*u(12))/dx		
  du(8) = (d70*u(1)+d71*u(2)+d72*u(3)+d73*u(4)+d74*u(5)+d75*u(6)+d76*u(7)+d77*u(8)+d78*u(9)+d79*u(10)+d710*u(11)+d711*u(12))/dx		
	
  do i=9,nx-8
  du(i) = (-dc4*u(i-4)-dc3*u(i-3)-dc2*u(i-2)-dc1*u(i-1) &
  &	   +dc4*u(i+4)+dc3*u(i+3)+dc2*u(i+2)+dc1*u(i+1))/dx
  end do		
	

  du(nx)   = -(d00*u(nx)+d01*u(nx-1)+d02*u(nx-2)+d03*u(nx-3)+d04*u(nx-4)+d05*u(nx-5)+d06*u(nx-6)+d07*u(nx-7)+d08*u(nx-8)+d09*u(nx-9)+d010*u(nx-10)+d011*u(nx-11))/dx
  du(nx-1) = -(d10*u(nx)+d11*u(nx-1)+d12*u(nx-2)+d13*u(nx-3)+d14*u(nx-4)+d15*u(nx-5)+d16*u(nx-6)+d17*u(nx-7)+d18*u(nx-8)+d19*u(nx-9)+d110*u(nx-10)+d111*u(nx-11))/dx
  du(nx-2) = -(d20*u(nx)+d21*u(nx-1)+d22*u(nx-2)+d23*u(nx-3)+d24*u(nx-4)+d25*u(nx-5)+d26*u(nx-6)+d27*u(nx-7)+d28*u(nx-8)+d29*u(nx-9)+d210*u(nx-10)+d211*u(nx-11))/dx
  du(nx-3) = -(d30*u(nx)+d31*u(nx-1)+d32*u(nx-2)+d33*u(nx-3)+d34*u(nx-4)+d35*u(nx-5)+d36*u(nx-6)+d37*u(nx-7)+d38*u(nx-8)+d39*u(nx-9)+d310*u(nx-10)+d311*u(nx-11))/dx		
  du(nx-4) = -(d40*u(nx)+d41*u(nx-1)+d42*u(nx-2)+d43*u(nx-3)+d44*u(nx-4)+d45*u(nx-5)+d46*u(nx-6)+d47*u(nx-7)+d48*u(nx-8)+d49*u(nx-9)+d410*u(nx-10)+d411*u(nx-11))/dx		
  du(nx-5) = -(d50*u(nx)+d51*u(nx-1)+d52*u(nx-2)+d53*u(nx-3)+d54*u(nx-4)+d55*u(nx-5)+d56*u(nx-6)+d57*u(nx-7)+d58*u(nx-8)+d59*u(nx-9)+d510*u(nx-10)+d511*u(nx-11))/dx		
  du(nx-6) = -(d60*u(nx)+d61*u(nx-1)+d62*u(nx-2)+d63*u(nx-3)+d64*u(nx-4)+d65*u(nx-5)+d66*u(nx-6)+d67*u(nx-7)+d68*u(nx-8)+d69*u(nx-9)+d610*u(nx-10)+d611*u(nx-11))/dx
  du(nx-7) = -(d70*u(nx)+d71*u(nx-1)+d72*u(nx-2)+d73*u(nx-3)+d74*u(nx-4)+d75*u(nx-5)+d76*u(nx-6)+d77*u(nx-7)+d78*u(nx-8)+d79*u(nx-9)+d710*u(nx-10)+d711*u(nx-11))/dx
  		
	end subroutine derivs48_b		
	END MODULE M_DERIVS
