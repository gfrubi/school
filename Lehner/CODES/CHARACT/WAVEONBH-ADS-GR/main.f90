!!! This code integrates the Einstein-scalar system in spherical
!!! symmetry, in ingoing Eddington-Finkelstein coordinates. See the
!!! LaTeX notes in the same directory for further description.

! Initial data consists of the scalar field, and boundary data is the
! (conserved) ADM mass at the initial time.

! All the gravitational variables, plus the time derivative of the
! scalar field, are integrated along constant-v (time) slices. Once
! those are obtained, the time integration is performed for the scalar
! field itself.

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
  integer :: i,j, Nv, Nsteps, Nspatial, Ntimelike
  integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
  logical :: ltrace	
  real*8 :: coeff(3), p0, dtp0, d2tp0,tangh, Axx,xm

  !set for debug or not
  ltrace = .true.

  !read pars
  call readpars
  !some for the excision stuff
  IMAX = N
  IMAXOLD = -1

  ! U will stand for the variables to study
  ! There are 7 variables total, 6 of which are integrated along
  ! spatial slices, 1 in time.
  Nv = 7
  Nspatial = 6
  Ntimelike = 1

  !allocate needed memmory
  !Unew will have the updated values of the fields
  !uold will have the old values of the fields
  !x, y are coords 

  !assume same extent for both coordinates
  allocate(Unew(Nv,N),Uold(Nv,N),x(N),tempout(N))


  !define coords and related. For simplicity
  !we do a straight uniform grid. from rhoin to rhout
  !but we label with i=1 the outer boundary point
  !we call this coordinate x
  ! x(1) = rhoin, x(N) = rhoout, so dx > 0

  dx= (rhoout-rhoin)/(N-1)
  dt = cfl * abs(dx)

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

  !ret = gft_out_full('phi',time, n, 'x', 1, x, Uold(iphi,:))
  !ret = gft_out_full('pi',time, n, 'x', 1, x, Uold(ipi,:))

  
  
  !print *, Uold(iphi,:)
  
  
  if(ltrace) print*, 'about to integrate'

  do ITE=1, nt
     !integrate
     
     !call evolution routine
     call evolve(unew,uold,dx,dt,time,x,Nv,Nspatial,N,lambda)

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
        ret = gft_out_full('pi',time, n, 'x', 1, x, Uold(ipi,:))
        ret = gft_out_full('s',time, n, 'x', 1, x, Uold(is,:))
        ret = gft_out_full('alpha',time, n, 'x', 1, x, Uold(ialpha,:))
        ret = gft_out_full('sigma',time, n, 'x', 1, x, Uold(isigma,:))
        ret = gft_out_full('beta',time, n, 'x', 1, x, Uold(ibeta,:))
        ret = gft_out_full('tau',time, n, 'x', 1, x, Uold(itau,:))
     end if


  end do

end program main
