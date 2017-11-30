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
  use M_apparent_horizon
  use M_resume

  implicit none

  !U will be vars
  real*8, dimension(:,:), allocatable :: Unew, Uold, aux_output_vec	
  real*8, dimension(:), allocatable :: x, tempout, dphir
  real*8 :: dx, dt, time
  integer :: i,j, Nv, Nc, Nsteps, Nspatial, Ntimelike, Nhor
  integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
  logical :: ltrace	
  real*8 :: coeff(3), p0, dtp0, d2tp0,tangh, Axx,xm
  real*8, dimension(:), allocatable :: rho_ah, rho_out_new, area_ah, darea_ah, horizon_out_vec
  integer :: N_total, ll
  real*8 :: chargeAss, chargeSFss, chargeAtt, chargeSFtt,chargeAts, chargeSFts, chargeSFur, chargeAur, charge_jr
  real*8 :: Pfluxv_int
  real*8, parameter :: PI = 3.1415926535897932



  !set for debug or not
  ltrace = .true.

  !read pars
  call readpars
  !some for the excision stuff
  IMAX = N
  IMAXOLD = N

  ! U will stand for the variables to study
  ! There are 7 variables total, 6 of which are integrated along
  ! spatial slices, 1 in time.
  Nv = 11
  Nspatial = 9
  Ntimelike = 2
  ! Number of auxiliary outputs
  Nc = 17
  ! Number of horizon quantities to output
  Nhor = 4


  !allocate needed memmory
  !Unew will have the updated values of the fields
  !uold will have the old values of the fields
  !x, y are coords 

  !assume same extent for both coordinates
  allocate(Unew(Nv,N),Uold(Nv,N),x(N),tempout(N),aux_output_vec(Nc,N), rho_ah(Nt), &
       &     rho_out_new(Nt), area_ah(Nt),darea_ah(Nt), dphir(N), &
       & horizon_out_vec(Nhor))


  !define coords and related. For simplicity
  !we do a straight uniform grid. from rhoin to rhout
  !but we label with i=1 the outer boundary point
  !we call this coordinate x
  ! x(1) = rhoin, x(N) = rhoout, so dx > 0

  print*, 'Opening files: status replace.'
  open(unit=10, file='area_time.dat',status='replace')
  open(unit=11, file='darea_time.dat',status='replace')
  open(unit=12, file='horizon_output.dat',status='replace')
  open(unit=111, file='dxphir_bdry.dat',status='replace')
  open(unit=112, file='stress-energy-int.dat',status='replace')
  open(unit=113, file='charge_jr.dat',status='replace')



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

  chargeAss = 0.0
  chargeSFss = 0.0
  chargeAtt = 0.0
  chargeSFtt = 0.0
  chargeAts= 0.0
  chargeSFts=0.0
  chargeSFur=0.0 
  chargeAur=0.0
  charge_jr=0.0
  Pfluxv_int=0.0



  if (resumelevel.eq.0) then
     call initial(Uold,dx,x,Nv,N,dt,lambda)
     time = initime
  else
     call resume(Uold,time,N)
  end if

  Unew = Uold

  if(ltrace) print*, 'about to integrate'

  !zero stuff
  horizon_out_vec = 0.0d0
  aux_output_vec = 0.0d0

  do ITE=1, nt
     !integrate

     !call evolution routine
     call evolve(unew,uold,aux_output_vec,dx,dt,time,x,Nv,Nspatial,N,lambda)


     if (ITE.eq.1.or.mod(ITE,freqhor).eq.0) then
        IMAXOLD = IMAX
        call apparent_horizon(unew(is,:),unew(isigma,:)&
             &,x,N_total,rho_ah(ITE),rho_out_new(ITE), area_ah(ITE))

        !need to fill values that might not exist if we went inwards
        if(IMAX.gt.IMAXOLD) then
           print*,'extrapol!', IMAX, IMAXOLD
           do i=IMAX,IMAXOLD
              unew(:,i) = 2.*unew(:,i-1) - unew(:,i-2)
           end do
        end if

        write(10,'(4(E14.7))') time, rho_ah(ITE), area_ah(ITE)/area_ah(1),0.25*sqrt(area_ah(ITE)/Pi)

     end if

     ! update the time
     !	dt = cfl * dx/(1.+abs(dtp0))
     time = time + dt

     ! Initial output of horizon quantities (after 1 time-step)
     if (ITE.eq.1.or.mod(ITE,freqhor).eq.0) then
        call horizon_output(Uold,horizon_out_vec,rho_ah(ITE),Nhor,N,Nv,x,dx)
        ! Columns are j_s, j_s, Tphi_ts, TA_ts
        write(12,'(5(E16.7))'), time, horizon_out_vec(ihorjs), &
             & horizon_out_vec(ihorTphits), horizon_out_vec(ihorTAts),horizon_out_vec(ihorcharge)
     end if


     !call screen output every now and then
     if (mod(ITE,freq).eq.0) then
        !          print*, 'call out', time

        write(*,'(4(E14.7))'), time,rho_ah(ITE), area_ah(ITE),0.25*sqrt(area_ah(ITE)/Pi)   

     end if

     if (mod(ITE,freqfreq).eq.0) then
        call derivs(unew(iphir,:),dphir(:),dx, N,1,IMAX)
	write(111,*), time, dphir(1)
     end if

     !some output 
     if (mod(ITE,freqsdf).eq.0) then
        ret = gft_out_full(trim(oprefix)//'phir',time, n, 'x', 1, x, Unew(iphir,:)*mask)
        ret = gft_out_full(trim(oprefix)//'phii',time, n, 'x', 1, x, Unew(iphii,:)*mask)
        ret = gft_out_full(trim(oprefix)//'pir',time, n, 'x', 1, x, Unew(ipir,:)*mask)
        ret = gft_out_full(trim(oprefix)//'pii',time, n, 'x', 1, x, Unew(ipii,:)*mask)
        ret = gft_out_full(trim(oprefix)//'s',time, n, 'x', 1, x, Unew(is,:)*mask)
        ret = gft_out_full(trim(oprefix)//'alpha',time, n, 'x', 1, x, Unew(ialpha,:)*mask)
        ret = gft_out_full(trim(oprefix)//'sigma',time, n, 'x', 1, x, Unew(isigma,:)*mask)
        ret = gft_out_full(trim(oprefix)//'beta',time, n, 'x', 1, x, Unew(ibeta,:)*mask)
        ret = gft_out_full(trim(oprefix)//'tau',time, n, 'x', 1, x, Unew(itau,:)*mask)
        ret = gft_out_full(trim(oprefix)//'w',time, n, 'x', 1, x, Unew(iw,:)*mask)
        ret = gft_out_full(trim(oprefix)//'z',time, n, 'x', 1, x, Unew(iz,:)*mask)
        ret = gft_out_full(trim(oprefix)//'mask',time, n, 'x', 1, x, mask)
        call aux_output(Uold,Unew,aux_output_vec,dx,dt,Nv,N,x)
        ret = gft_out_full(trim(oprefix)//'jn',time, n, 'x', 1, x, aux_output_vec(ijn,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'js',time, n, 'x', 1, x, aux_output_vec(ijs,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'Tphitt',time, n, 'x', 1, x, aux_output_vec(iTphitt,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'Tphits',time, n, 'x', 1, x, aux_output_vec(iTphits,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'Tphiss',time, n, 'x', 1, x, aux_output_vec(iTphiss,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'TAtt',time, n, 'x', 1, x, aux_output_vec(iTAtt,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'TAts',time, n, 'x', 1, x, aux_output_vec(iTAts,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'TAss',time, n, 'x', 1, x, aux_output_vec(iTAss,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'kretschmann',time, n, 'x', 1, x, aux_output_vec(ikretsch,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'rho2d+Sigma',time, n, 'x', 1, x, aux_output_vec(idplussigma,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'res1',time, n, 'x', 1, x, aux_output_vec(ires1,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'res2',time, n, 'x', 1, x, aux_output_vec(ires2,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'charge',time, n, 'x', 1, x, aux_output_vec(icharge,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'misner',time, n, 'x', 1, x, aux_output_vec(imisner,:)*maskhor)
        ret = gft_out_full(trim(oprefix)//'Pfluxv',time, n, 'x', 1, x, aux_output_vec(iPfluxv,:)*maskhor)
        call horizon_output(Uold,horizon_out_vec,rho_ah(ITE),Nhor,N,Nv,x,dx)
        ! Columns are j_s, j_s, Tphi_ts, TA_ts
        write(12,'(5(E16.7))'), time, horizon_out_vec(ihorjs), &
             & horizon_out_vec(ihorTphits), horizon_out_vec(ihorTAts), horizon_out_vec(ihorcharge)

        !integrate, brute-forcely now the 'charges'
        do ll = 10, n-1
           chargeAss = chargeAss +   dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTAss,ll)    *maskhor(ll)
           chargeSFss = chargeSFss + dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTphiss,ll)  *maskhor(ll)
           chargeAtt = chargeAtt +   dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTAtt,ll)    *maskhor(ll)
           chargeSFtt = chargeSFtt + dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTphitt,ll)  *maskhor(ll)
           chargeAts = chargeAts +   dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTAts,ll)    *maskhor(ll)
           chargeSFts = chargeSFts + dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTphits,ll)  *maskhor(ll)
           chargeAur = chargeAur +   dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTAvr,ll)    *maskhor(ll)
           chargeSFur = chargeSFur + dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iTphivr,ll)  *maskhor(ll)
           charge_jr = charge_jr - dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(ijs,ll)  *maskhor(ll)
           Pfluxv_int = Pfluxv_int - dx/x(ll)**2 * (unew(isigma,ll) + 1./x(ll))**2* aux_output_vec(iPfluxv,ll)  *maskhor(ll)
	end do
        write(112,'(10(E16.7))') time,chargeAtt,chargeAts,chargeAss,chargeAur,chargeSFtt,&
             & chargeSFts,chargeSFss,chargeSFur,Pfluxv_int
        write(113,'(2(E16.7))') time,charge_jr

	chargeAss = 0.0d0
        chargeSFss = 0.0d0
	chargeAtt = 0.0d0
        chargeSFtt = 0.0d0
	chargeAts = 0.0d0
        chargeSFts = 0.0d0
	chargeAur = 0.0d0
        chargeSFur = 0.0d0
        charge_jr = 0.0d0
        Pfluxv_int = 0.0d0
     end if

     ! Save the new values in the old array
     uold = unew

  end do

  call area_evolution(area_ah,dt,darea_ah)
  write(11,'(2(E14.7))') (j*dt,darea_ah(j+1), j = 0, Nt-1)
  close(10)
  close(11)
  close(12)
  close(111)
  close(112)
  close(113)

end program main
