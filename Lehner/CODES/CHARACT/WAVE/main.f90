!this code will integrate the wave equation 2 phi,ur + phi,rr =0
!using it in first order form. ie. we define g = phi,r
!and so the equations of motion become
!phi,r = g and g,u = -g,r/2
!for boundary conditions we will chose no incoming modes at the 
!right boundary r (or x as it will be called) = rhoout. Since the incoming
!mode corresponds to phi, this means we will choose phi=0
!in this code we will also need an inner boundary condition
!and also for simplicity we will adopt no incoming there, 
!which means g=0 at the inner boundary at rhoin


	program main
	use M_initial
	use M_evolve
	use M_pars
	use M_derivs
	
	implicit none
	
!U will be vars
	real*8, dimension(:,:), allocatable :: Unew, Uold	
	real*8, dimension(:), allocatable :: x, tempout
	real*8 :: dx, dt, time
	integer :: i,j, Nv, Nsteps
        integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
	logical :: ltrace	
        real*8 :: coeff(3), p0, dtp0, d2tp0, Axx,xm

!set for debug or not
	ltrace = .true.

!read pars
        call readpars
!some for the excision stuff
        IMAX = N
        IMAXOLD = -1

!U will stand for the variables to study
	Nv = 6

!allocate needed memory
!Unew will have the updated values of the fields
!uold will have the old values of the fields
!x are coords 

!assume same extent for both coordinates
	allocate(Unew(Nv,N),Uold(Nv,N),x(N),tempout(N))
	
		
!define coords and related. For simplicity
!we do a straight uniform grid. from rhoin to rhout
!but we label with i=1 the outer boundary point
!we call this coordinate x 
	
	dx= (rhoout-rhoin)/(N-1)
	dt = cfl * dx
	
!straight out
	do i = 1, N
	x(i) =  rhoin + (i-1)*dx
	end do
	
!INITIALIZE
!pars for talking from the boundary
	a2i = 0.0
        p2 = 0.0
	a2iold = -1.0
        p2old = 0.0

	call initial(Uold,dx,x,Nv,N,dt,lambda)
        time = initime	
	Unew = Uold

	if(ltrace) print*, 'about to integrate'
	
	   do ITE=1, nt
	    !integrate

!call evolution routine
	    call evolve(unew,uold,dx,dt,time,x,Nv,N,lambda)

           !update the time and save the new values
           !in the old array.

!this was here adjust time step, to react to rapid changes
!but it is not needed here so commenting it out

!	dt = cfl * dx/(1.+abs(dtp0))
            time = time + dt

            uold = unew

!call screen output every now and then
          if (mod(ITE,freq).eq.0) then
!          print*, 'call out', time
           print*, time, p2, a2
 	  end if

!some output 
          if (mod(ITE,freqsdf).eq.0) then
          ret = gft_out_full('phi',time, n, 'x', 1, x, Uold(iphi,:))
          ret = gft_out_full('gg',time, n, 'x', 1, x, Uold(ig,:))

          end if


	   end do

	end
