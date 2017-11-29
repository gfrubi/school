	MODULE M_PARS
	implicit none
	
!define parameters !!!!!!!!!!!!!!!!!!!!!!!
	real*8 ::  rhoin, rhoout, cfl, disip, gamma, amp,lambda,initime
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
	&		gamma, amp,lambda, initime, rhoin, rhoout

!read params !!!!!!!!!!!!!!!!!!!!!!!


   open (unit = 10, file = "pars.in", status = "old" )
   read (unit = 10, nml = pars_input)
   close(unit = 10)


!we add also some pars to label fields so that
!we can later call the fields without having to remember
!what index corresponds to what field
!since we only use the first two the others do not matter
!here using iphi to 1 giving the field phi
!and ig the point 2 giving the field g
	iphi = 1
        ig = 2
        ipi = 3
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
