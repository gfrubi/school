MODULE M_PARS
  implicit none

  !define parameters !!!!!!!!!!!!!!!!!!!!!!!
  real*8 ::  rhoin, rhoout, cfl, disip, gamma, amp,lambda,initime,mbh, rads, rh
  real*8 :: rhoh, rhoads
  integer :: N, Nt, freqsdf, freq, freqhor, bc, derorder,dissorder,STMETH

  integer :: iphi, ig,ipi, ib,ibeta,ia, ialp, IMAX, IMAXOLD, ITE
  real*8 :: a2i, p2, p2old, a2iold, a2, area_rad, null_rad

  real*8, dimension(:), allocatable :: mask

CONTAINS

  !these are some generic parameters to use when the time calls for them	
  subroutine readpars
    implicit none
    namelist /pars_input/ N, Nt, freqsdf, freq,freqhor, cfl, disip, &
         &		derorder,dissorder, STMETH, &
         &		gamma, amp,lambda, initime, rhoin, rhoout, rhoh, rhoads

    !read params !!!!!!!!!!!!!!!!!!!!!!!


    open (unit = 10, file = "pars.in", status = "old" )
    read (unit = 10, nml = pars_input)
    close(unit = 10)

    ! Black hole horizon
!!$    rh = (-(3.**0.6666666666666666*Lads**2.) + 3.**0.3333333333333333* &
!!$         &     (Lads**2.*(9.*mbh + Sqrt(3.)*Sqrt(Lads**2. + 27.*mbh**2.)))**0.6666666666666666)/ &
!!$         &  (3.*(Lads**2.*(9.*mbh + Sqrt(3.)*Sqrt(Lads**2. + 27.*mbh**2.)))**0.3333333333333333)

    ! Radius conversions
!!$    rhoh = 1./rh
    rh = 1./rhoh
    rads = 1./rhoads

    ! Black hole mass
    mbh = 1/(2.*rhoh) + 1/(2.*(rhoh**3.)*(rads**2.))
    
    ! Output some parameters
    print *, "Black hole mass:",mbh
    print *, "AdS radius, r = ",rads, " rho =", rhoads
    print *, "Horizon radius, r =", rh, " rho =", rhoh
    print *, "Domain, rhoin =", rhoin, " to rhoout =", rhoout

    !we add also some pars to label fields so that
    !we can later call the fields without having to remember
    !what index corresponds to what field
    !since we only use the first two the others do not matter
    !here using iphi to 1 giving the field phi
    !and ig the point 2 giving the field g
    iphi = 1
    ipi = 2
    ig = 3
    ibeta = 4
    ia = 5
    ialp = 6

    IMAX = N

    !mask is a convenient (at some point) array that may be needed to mask
    !things eventually inside the horizon if one
    !wanted to keep track of things
    allocate(mask(N))

    mask = 1.
    a2 = -1.
    a2i = 0.
    a2iold = 0.
    p2 = 0.
    p2old = 0.

  end subroutine readpars


end module M_pars
