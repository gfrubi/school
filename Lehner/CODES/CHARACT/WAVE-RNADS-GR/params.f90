MODULE M_PARS
  implicit none

  !define parameters !!!!!!!!!!!!!!!!!!!!!!!
  real*8 ::  rhoin, rhoout, cfl, disip, gamma, amp,lambda,initime,madm,qadm, rads, rh, q
  real*8 :: rhoh, rhoads
  integer :: N, Nt, freqsdf, freq, freqhor, freqfreq, bc, derorder,dissorder,STMETH,resumelevel

  integer :: isigma, itau, is, ialpha, ibeta, iw, iz, ipir, ipii, iphir, iphii, &
     &       ijn, ijs, iTphitt, iTphits, iTAtt, iTAts, icharge, iTAss,iTphiss, &
     &       imisner,ikretsch,idplussigma, iTAvr, iTphivr, iPfluxv, &
     &       ires1, ires2, ihorjs, ihorTphits, ihorTAts, ihorcharge, IMAX, IMAXOLD, ITE
  real*8 :: a2i, p2, p2old, a2iold, a2, area_rad, null_rad
  logical :: scalaronly
  character(len=32) :: iprefix, oprefix


  real*8, dimension(:), allocatable :: mask, maskhor

CONTAINS

  !these are some generic parameters to use when the time calls for them	
  subroutine readpars
    implicit none
    namelist /pars_input/ N, Nt, freqsdf, freq,freqhor, freqfreq, cfl, disip, &
         &		derorder,dissorder, STMETH, madm, qadm, &
         &		gamma, amp,lambda, initime, rhoin, rhoout, rhoads, scalaronly,q, &
         &              resumelevel, iprefix, oprefix

    !read params !!!!!!!!!!!!!!!!!!!!!!!


    ! Default values
    resumelevel = 0 ! Start with initial data instead of resuming
    iprefix = "" ! Input file prefix
    oprefix = "" ! Output file prefix

    ! Load from pars.in file
    open (unit = 10, file = "pars.in", status = "old" )
    read (unit = 10, nml = pars_input)
    close(unit = 10)



    ! AdS radius conversion
    rads = 1./rhoads

    ! Black hole horizon (approximate; ignores scalar field)
    rh = (Sqrt(3.)*Sqrt(-2.*rads**2. +  &
         & (3.*qadm**2.*rads + rads**3.)/ &
         & (54.*madm**2.*rads - 9.*qadm**2.*rads + rads**3. +  &
         & 3.*Sqrt(3.)*Sqrt(-qadm**6. +  &
         & 2.*(54.*madm**4. - 18.*madm**2.*qadm**2. + qadm**4.)*rads**2. +  &
         & (4.*madm**2. - qadm**2.)*rads**4.))**0.3333333333333333 +  &
         & (9.*(6.*madm**2. - qadm**2.)*rads**4. + rads**6. +  &
         & 3.*Sqrt(3.)*rads**3.*Sqrt(-qadm**6. +  &
         & 2.*(54.*madm**4. - 18.*madm**2.*qadm**2. + qadm**4.)*rads**2. +  &
         & (4.*madm**2. - qadm**2.)*rads**4.))**0.3333333333333333) +  &
         & 3.*Sqrt((-4.*rads**2.)/3. -  &
         & (rads*(3.*qadm**2. + rads**2.))/ &
         & (3.*(54.*madm**2.*rads - 9.*qadm**2.*rads + rads**3. +  &
         & 3.*Sqrt(3.)*Sqrt(-qadm**6. +  &
         & 2.*(54.*madm**4. - 18.*madm**2.*qadm**2. + qadm**4.)*rads**2. +  &
         & (4.*madm**2. - qadm**2.)*rads**4.))**0.3333333333333333) -  &
         & (rads*(54.*madm**2.*rads - 9.*qadm**2.*rads + rads**3. +  &
         & 3.*Sqrt(3.)*Sqrt(-qadm**6. +  &
         & 2.*(54.*madm**4. - 18.*madm**2.*qadm**2. + qadm**4.)*rads**2. +  &
         & (4.*madm**2. - qadm**2.)*rads**4.))**0.3333333333333333)/3. +  &
         & (4.*Sqrt(3.)*madm*rads**2.)/ &
         & Sqrt(-2.*rads**2. + (3.*qadm**2.*rads + rads**3.)/ &
         & (54.*madm**2.*rads - 9.*qadm**2.*rads + rads**3. +  &
         & 3.*Sqrt(3.)*Sqrt(-qadm**6. +  &
         & 2.*(54.*madm**4. - 18.*madm**2.*qadm**2. + qadm**4.)*rads**2. +  &
         & (4.*madm**2. - qadm**2.)*rads**4.))**0.3333333333333333 +  &
         & (9.*(6.*madm**2. - qadm**2.)*rads**4. + rads**6. +  &
         & 3.*Sqrt(3.)*rads**3.* &
         & Sqrt(-qadm**6. +  &
         & 2.*(54.*madm**4. - 18.*madm**2.*qadm**2. + qadm**4.)*rads**2. +  &
         & (4.*madm**2. - qadm**2.)*rads**4.))**0.3333333333333333)))/6.
    rhoh = 1./rh
      
    ! Output some parameters
    print *, "charge parameter:",q
    print *, "ADM mass:",madm
    print *, "ADM charge:",qadm
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
    iw = 6
    iz = 7
    ipir = 8
    ipii = 9

    ! Fields that evolve in time only
    iphir = 10
    iphii = 11

    ! Auxiliary output components
    ijn = 1
    ijs = 2
    ikretsch = 3
    idplussigma = 4
    ires1 = 5
    ires2 = 6
    iTphitt = 7
    iTphits = 8
    iTAtt = 9
    iTAts = 10
    icharge = 11
    iTphiss = 12
    iTAss = 13
    imisner = 14
    iTAvr = 15
    iTphivr = 16
    iPfluxv = 17

    ! Horizon output components
    ihorjs = 1
    ihorTAts = 2
    ihorTphits = 3
    ihorcharge = 4

    IMAX = N

    !mask is a convenient (at some point) array that may be needed to mask
    !things eventually inside the horizon if one
    !wanted to keep track of things
    allocate(mask(N),maskhor(N))

    mask = 1.
    maskhor = 1.
    a2 = -1.
    a2i = 0.
    a2iold = 0.
    p2 = 0.
    p2old = 0.

  end subroutine readpars


end module M_pars
