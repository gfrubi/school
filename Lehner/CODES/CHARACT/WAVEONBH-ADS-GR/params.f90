MODULE M_PARS
  implicit none

  !define parameters !!!!!!!!!!!!!!!!!!!!!!!
  real*8 ::  rhoin, rhoout, cfl, disip, gamma, amp,lambda,initime,madm, rads, rh
  real*8 :: rhoh, rhoads
  integer :: N, Nt, freqsdf, freq, freqhor, bc, derorder,dissorder,STMETH

  integer :: isigma, itau, is, ialpha, ibeta, ipi, iphi, IMAX, IMAXOLD, ITE
  real*8 :: a2i, p2, p2old, a2iold, a2, area_rad, null_rad
  logical :: scalaronly

  real*8, dimension(:), allocatable :: mask

CONTAINS

  !these are some generic parameters to use when the time calls for them	
  subroutine readpars
    implicit none
    namelist /pars_input/ N, Nt, freqsdf, freq,freqhor, cfl, disip, &
         &		derorder,dissorder, STMETH, madm, &
         &		gamma, amp,lambda, initime, rhoin, rhoout, rhoads, scalaronly

    !read params !!!!!!!!!!!!!!!!!!!!!!!


    open (unit = 10, file = "pars.in", status = "old" )
    read (unit = 10, nml = pars_input)
    close(unit = 10)



    ! AdS radius conversion
    rads = 1./rhoads

    ! Black hole horizon (approximate; ignores scalar field)
    rh = (-(3.**0.6666666666666666*rads**2.) + 3.**0.3333333333333333* &
         &     (rads**2.*(9.*madm + Sqrt(3.)*Sqrt(rads**2. + 27.*madm**2.)))**0.6666666666666666)/ &
         &  (3.*(rads**2.*(9.*madm + Sqrt(3.)*Sqrt(rads**2. + 27.*madm**2.)))**0.3333333333333333)
    rhoh = 1./rh
      
    ! Output some parameters
    print *, "ADM mass:",madm
    print *, "AdS radius, r = ",rads, " rho =", rhoads
    print *, "Approximate horizon radius (neglecting scalar), r =", rh, " rho =", rhoh
    print *, "Domain, rhoin =", rhoin, " to rhoout =", rhoout
    print *, "Scalar field only (no GR):", scalaronly

    !we add also some pars to label fields so that
    !we can later call the fields without having to remember
    !what index corresponds to what field

    ! Fields that are evolved by spatial integration
    isigma = 1
    itau = 2
    is = 3
    ialpha = 4
    ibeta = 5
    ipi = 6

    ! Fields that evolve in time only
    iphi = 7
    
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
