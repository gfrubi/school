!this file contains the RHSs of the equations.
!divided in those needed for the time integration
!and those for the space integration


! Useful regexp replace... \([^\.][[:digit:]]+\)\([^\.]\) -> \1.\2

MODULE M_RHS
  implicit none

contains


  subroutine rhs_time(u,du,source,dx,N,x,time)
    use m_pars, only : amp,IMAX, isigma, itau, is, ialpha, ibeta, iw, iz, ipir, ipii, iphir, iphii, &
         & lambda, p2, a2, madm, qadm, rhoads, rads, scalaronly
    implicit none
    real*8, dimension(:,:), intent(in) :: u, du
    real*8, dimension(:,:), intent(inout) :: source
    real*8, dimension(:), intent(in) :: x
    real*8 :: dx, v,time, uphir,uphii,usigma,utau,us,ualpha,ubeta,uw,uz,upir,upii,dxuphir,dxuphii
    integer :: i, N
    logical :: ltrace
    real*8 :: tangh,p0,dt_p0,d2t_p0,d3t_p0, xm, Axx, aa

    !for debug
    ltrace=.false.

    !here in case we want to keep messy things in 
    !a separate module	
    source(:,:) = 0.0

    ! Boundary evolution
    source(iphir,1) = 0.0d0
    source(iphii,1) = 0.0d0

    ! Evolution away from boundary
    do i=2, IMAX
       uphir = u(iphir,i)
       dxuphir = du(iphir,i)
       uphii = u(iphii,i)
       dxuphii = du(iphii,i)
       usigma = u(isigma,i)
       utau = u(itau,i)
       us = u(is,i)
       ualpha = u(ialpha,i)
       ubeta = u(ibeta,i)
       uw = u(iw,i)
       uz = u(iz,i)
       upir = u(ipir,i)
       upii = u(ipii,i)

       xm = x(i)

       source(iphir,i) = (2.*rads**2.*upir +  &
            & (2.*uphir + dxuphir*xm)* &
            & (1. + 2.*usigma*xm + rads**2.*xm**2. +  &
            & rads**2.*ualpha*xm**2. + usigma**2.*xm**2.))/ &
            & (2.*rads**2.*xm)

       source(iphii,i) = (2.*rads**2.*upii +  &
            & (2.*uphii + dxuphii*xm)* &
            & (1. + 2.*usigma*xm + rads**2.*xm**2. +  &
            & rads**2.*ualpha*xm**2. + usigma**2.*xm**2.))/ &
            & (2.*rads**2.*xm)
    end do
    !     source(iphi,1) = 0d0 !.0.5d0*u(ig,1)
    !	source(ig,IMAX) = -du(ig,IMAX)
    !	print*,'source',source(iphi,1)
    !	source(ig,IMAX) = -0.5d0*(du(iphi,IMAX)-du(iphi,IMAX-1))/dx


  end subroutine rhs_time

  subroutine rhs_space(uguess,u,du,d2u,source,dx,Nv,Nspatial,N,xm,time)
    use m_pars, only : amp, gamma,isigma, itau, is, ialpha, ibeta, iw, iz, ipir, ipii, iphir, iphii, lambda,a2,p2,rhoads,madm,qadm,rads,scalaronly,q
    implicit none
    real*8, dimension(:), intent(in) :: u, du
    real*8 :: time, xm

    real*8, dimension(Nv) :: source, d2u
    real*8, dimension(Nv) :: uguess
    real*8 :: dx, v
    integer :: i, j, Nv, Nspatial, N
    logical :: ltrace

    real*8, dimension(Nv) :: At,Arr,drArr,Ar,A0,J0
    real*8 :: tangh,p0,dt_p0,d2t_p0,d3t_p0
    real*8 :: uphir,uphii,usigma,utau,us,ualpha,ubeta,uw,uz,upir,upii,bet,bb,alp,aa
    real*8 :: dx_ppi,dxuphir,dxuphii,dx_bet,dx_bb,dx_alp,dx_aa
    real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,t0,rho
    real*8 :: RR(3), lnxm, mod
    !for debug
    ltrace=.false.


    !define first all that is needed, notice phi is evolved
    !already so they come from u, the rest from uguess

    uphir = u(iphir)
    dxuphir = du(iphir)

    uphii = u(iphii)
    dxuphii = du(iphii)    

    usigma = uguess(isigma)
    utau = uguess(itau)
    us = uguess(is)
    ualpha = uguess(ialpha)
    ubeta = uguess(ibeta)
    uw = uguess(iw)
    uz = uguess(iz)
    upir = uguess(ipir)
    upii = uguess(ipii)



    source = 0.

    ! First calculate source for scalar field

    ! We have a special case at rho = 0 end point
    if (xm.lt.(0.1*dx)) then
       source(ipir) = -1.5*dxuphir/rads**2.

       source(ipii) = -1.5*dxuphii/rads**2.
    else
       source(ipir) = -(2.*uphir + dxuphir*xm + 4.*uphir*usigma*xm +  &
            & 2.*rads**2.*uphir*xm**2. + 4.*rads**2.*uphir*us*xm**2. +  &
            & 2.*dxuphir*usigma*xm**2. + 2.*uphir*usigma**2.*xm**2. +  &
            & q*rads**2.*uphii*uz*xm**2. + dxuphir*rads**2.*xm**3. +  &
            & 2.*dxuphir*rads**2.*us*xm**3. +  &
            & dxuphir*usigma**2.*xm**3. +  &
            & q*rads**2.*uphii*usigma*uz*xm**3. +  &
            & 2.*rads**2.*upir*xm*(usigma + utau*xm) +  &
            & 2.*q*rads**2.*uw*xm* &
            & (dxuphii*xm*(1. + usigma*xm) +  &
            & uphii*(1. + 2.*usigma*xm + utau*xm**2.)))/ &
            & (2.*rads**2.*xm*(1. + usigma*xm))

       source(ipii) = (-2.*rads**2.*upii*xm*(usigma + utau*xm) -  &
            & 2.*uphii*(1. + 2.*usigma*xm + rads**2.*xm**2. +  &
            & 2.*rads**2.*us*xm**2. + usigma**2.*xm**2.) +  &
            & xm*(q*rads**2.*uphir*uz*xm*(1. + usigma*xm) -  &
            & dxuphii*(1. + 2.*usigma*xm + rads**2.*xm**2. +  &
            & 2.*rads**2.*us*xm**2. + usigma**2.*xm**2.) +  &
            & 2.*q*rads**2.*uw* &
            & (dxuphir*xm*(1. + usigma*xm) +  &
            & uphir*(1. + 2.*usigma*xm + utau*xm**2.))))/ &
            & (2.*rads**2.*xm*(1. + usigma*xm))
    end if

    ! To turn off coupling of the scalar, just set pi, phi equal to
    ! zero temporarily, before calculating metric and Maxwell fields
    if (scalaronly) then
       uphir = 0.
       uphii = 0.
       dxuphir = 0.
       dxuphii = 0.
       upii = 0.
       upir = 0.
    end if

    ! End point
    if (xm.lt.(0.1*dx)) then
       !print *, "end point treatment"
       source(is) = -madm
       source(ialpha) = -2.*madm
       source(ibeta) = 0.5*(qadm**2.)
       source(isigma) = 0.
       source(itau) = 0.
       source(iw) = qadm
       source(iz) = 0.
    else

       ! The s field receives a special implicit treatment. Zero out a
       ! contribution that goes as s/rho.
       mod=0.
       source(is) = (rads**2.*uz**2.*xm**2.*(1. + usigma*xm)**2. -  &
            & 4.*utau*(3. + 6.*usigma*xm + rads**2.*xm**2. +  &
            & 3.*usigma**2.*xm**2.) - 8.*rads**2.*us*(-1.*mod + utau*xm**2.) &
            & )/(8.*rads**2.*xm*(1. + usigma*xm))

       source(ialpha) = ubeta

       source(ibeta) = (2.*utau - 2.*rads**2.*ubeta*xm + 4.*usigma*utau*xm +  &
            & 4.*uphii**2.*xm**2. + 4.*uphir**2.*xm**2. +  &
            & 4.*rads**2.*uphii*upii*xm**2. +  &
            & 4.*rads**2.*uphir*upir*xm**2. -  &
            & 4.*rads**2.*ubeta*usigma*xm**2. - 2.*rads**2.*utau*xm**2. +  &
            & 2.*usigma**2.*utau*xm**2. - 2.*utau**2.*xm**2. +  &
            & 4.*dxuphii*uphii*xm**3. + 4.*dxuphir*uphir*xm**3. +  &
            & 2.*dxuphii*rads**2.*upii*xm**3. +  &
            & 2.*dxuphir*rads**2.*upir*xm**3. +  &
            & 16.*uphii**2.*usigma*xm**3. + 16.*uphir**2.*usigma*xm**3. +  &
            & 8.*rads**2.*uphii*upii*usigma*xm**3. +  &
            & 8.*rads**2.*uphir*upir*usigma*xm**3. -  &
            & 2.*rads**2.*ubeta*usigma**2.*xm**3. -  &
            & 4.*usigma*utau**2.*xm**3. + dxuphii**2.*xm**4. +  &
            & dxuphir**2.*xm**4. + 16.*dxuphii*uphii*usigma*xm**4. +  &
            & 16.*dxuphir*uphir*usigma*xm**4. +  &
            & 4.*dxuphii*rads**2.*upii*usigma*xm**4. +  &
            & 4.*dxuphir*rads**2.*upir*usigma*xm**4. +  &
            & 24.*uphii**2.*usigma**2.*xm**4. +  &
            & 24.*uphir**2.*usigma**2.*xm**4. +  &
            & 4.*rads**2.*uphii*upii*usigma**2.*xm**4. +  &
            & 4.*rads**2.*uphir*upir*usigma**2.*xm**4. -  &
            & 2.*usigma**2.*utau**2.*xm**4. +  &
            & 2.*dxuphir*q*rads**2.*uphii*uw*xm**4. -  &
            & 2.*dxuphii*q*rads**2.*uphir*uw*xm**4. +  &
            & 4.*dxuphii**2.*usigma*xm**5. +  &
            & 4.*dxuphir**2.*usigma*xm**5. +  &
            & 24.*dxuphii*uphii*usigma**2.*xm**5. +  &
            & 24.*dxuphir*uphir*usigma**2.*xm**5. +  &
            & 2.*dxuphii*rads**2.*upii*usigma**2.*xm**5. +  &
            & 2.*dxuphir*rads**2.*upir*usigma**2.*xm**5. +  &
            & 16.*uphii**2.*usigma**3.*xm**5. +  &
            & 16.*uphir**2.*usigma**3.*xm**5. +  &
            & 4.*dxuphir*q*rads**2.*uphii*usigma*uw*xm**5. -  &
            & 4.*dxuphii*q*rads**2.*uphir*usigma*uw*xm**5. +  &
            & 6.*dxuphii**2.*usigma**2.*xm**6. +  &
            & 6.*dxuphir**2.*usigma**2.*xm**6. +  &
            & 16.*dxuphii*uphii*usigma**3.*xm**6. +  &
            & 16.*dxuphir*uphir*usigma**3.*xm**6. +  &
            & 4.*uphii**2.*usigma**4.*xm**6. +  &
            & 4.*uphir**2.*usigma**4.*xm**6. +  &
            & 2.*dxuphir*q*rads**2.*uphii*usigma**2.*uw*xm**6. -  &
            & 2.*dxuphii*q*rads**2.*uphir*usigma**2.*uw*xm**6. +  &
            & 4.*dxuphii**2.*usigma**3.*xm**7. +  &
            & 4.*dxuphir**2.*usigma**3.*xm**7. +  &
            & 4.*dxuphii*uphii*usigma**4.*xm**7. +  &
            & 4.*dxuphir*uphir*usigma**4.*xm**7. +  &
            & dxuphii**2.*usigma**4.*xm**8. +  &
            & dxuphir**2.*usigma**4.*xm**8. +  &
            & rads**2.*uz**2.*xm**2.*(1. + usigma*xm)**2. -  &
            & 4.*rads**2.*us*(-1. + utau*xm**2.))/ &
            & (rads**2.*xm**2.*(1. + usigma*xm)**2.)

!!$       print *, "usigma", usigma
!!$       print *, "utau", utau
!!$       print *, "us", us
!!$       print *, "ualpha", ualpha
!!$       print *, "ubeta", ubeta
!!$       print *, "uw", uw
!!$       print *, "uz", uz
!!$       print *, "xm", xm
!!$       print *, "source(ibeta)", source(ibeta)

       source(isigma) = utau

       source(itau) = -(4.*utau + xm**2.*(1. + usigma*xm)* &
            & (4.*uphii**2. + 4.*uphir**2. + 4.*dxuphii*uphii*xm +  &
            & 4.*dxuphir*uphir*xm +  &
            & (dxuphii**2. + dxuphir**2.)*xm**2.))/(2.*xm)

       source(iw) = uz

       source(iz) = (-2.*(q*(-(dxuphir*uphii) + dxuphii*uphir)*xm**2.* &
            & (1. + usigma*xm) + uz*(usigma + utau*xm)))/ &
            & (1. + usigma*xm)

       ! Special correction for some of the fields near the boundary !
!!$       if (xm.lt.100.1*dx) then
!!$
!!$          source(ialpha) = -2.*madm + (qadm**2.)*xm/2.
!!$
!!$          source(ibeta) = (qadm**2.)/2.
!!$
!!$          source(is) = -madm + (qadm**2.)*xm/4.
!!$
!!$       end if
    
    end if

  end subroutine rhs_space


  !
  ! Return the rhs of the s equation, not including the term
  ! proportional to s/rho. This is needed for the special implicit
  ! integration of s.
  !
  function rhs_s(uguess,xm) result(s_result)
    use m_pars, only: rads, isigma, itau, is, ialpha, ibeta, iw, iz
    implicit none
    real*8, dimension(:), intent(in) :: uguess
    real*8, intent(in) :: xm
    real*8 :: s_result
    real*8 :: mod,usigma,utau,us,ualpha,ubeta,uw,uz

    usigma = uguess(isigma)
    utau = uguess(itau)
    us = uguess(is)
    ualpha = uguess(ialpha)
    ubeta = uguess(ibeta)
    uw = uguess(iw)
    uz = uguess(iz)
    
    !  Zero out a contribution that goes as s/rho.
    mod=0.
    s_result = (rads**2.*uz**2.*xm**2.*(1. + usigma*xm)**2. -  &
         & 4.*utau*(3. + 6.*usigma*xm + rads**2.*xm**2. +  &
         & 3.*usigma**2.*xm**2.) - 8.*rads**2.*us*(-1.*mod + utau*xm**2.) &
         & )/(8.*rads**2.*xm*(1. + usigma*xm))
    
  end function rhs_s
  

  !
  ! Calculate stess-energy, charge-current, Kretchmann, \rho^2
  ! d_+\Sigma (apparent horizon where this vanishes), and two
  ! independent residuals (unused equations of motion).
  !
  subroutine aux_output_rhs(u,dxu,dvu,output_vec,dx,xm)
    use m_pars, only : isigma, itau, is, ialpha, ibeta, iw, iz, ipir, &
         & ipii, iphii, iphir, ijn, ijs, iTphitt, iTphits, iTAtt, iTAts, ikretsch, iTAvr, iTphivr, &
         & idplussigma, ires1, ires2, rads, scalaronly, q, icharge,N, iTAss, iTphiss, imisner, &
         & iPfluxv
    implicit none
    real*8, dimension(:), intent(in) :: u, dxu, dvu
    real*8, dimension(:), intent(inout) :: output_vec
    real*8, intent(in) :: xm,dx
    real*8 :: uphir,uphii,usigma,utau,us,ualpha,ubeta,uw,uz,upir,upii
    real*8 :: dxuphir,dxuphii,dxutau,dxus,dxubeta,dxuz
    real*8 :: dvusigma,dvus,dvuz
    real*8 :: Tphi_vr, Tphi_vv, Tphi_rr, TA_vr, TA_vv, TA_rr, J_v, J_r
    real*8 :: t_v, t_r, n_v, n_r, s_v, s_r, u_v, u_r, xm_reg, DR, dplus

!we use A, Sigma too many times, let us define some temp array for those
    real*8 :: Afield, Sigmafield


    ! Needed variables
    uphir = u(iphir)
    dxuphir = dxu(iphir)

    uphii = u(iphii)
    dxuphii = dxu(iphii)    

    usigma = u(isigma)
    dvusigma = dvu(isigma)

    utau = u(itau)
    dxutau = dxu(itau)

    us = u(is)
    dxus = dxu(is)
    dvus = dvu(is)

    ualpha = u(ialpha)

    ubeta = u(ibeta)
    dxubeta = dxu(ibeta)

    uw = u(iw)

    uz = u(iz)
    dxuz = dxu(iz)
    dvuz = dvu(iz)

    upir = u(ipir)
    upii = u(ipii)


    ! Define vectors (up-index):

    ! Use a regulator, since A blows up at infinity (ad hoc). If
    ! fields go smoothly to 0 in the end, this is fine.
    if (xm.lt.0.0001*dx) then
       xm_reg = 0.0001*dx
    else
       xm_reg = xm
    end if


    Afield = (ualpha + (usigma + 1./xm_reg)**2./rads**2. + 1.)
    Sigmafield = usigma + 1./xm_reg
    dplus = us + 0.5d0 + usigma**2*0.5d0/rads**2

    ! t^a = (d/dv)^a + (A-1)/2 (d/dr)^a
    t_v = 1.
    t_r = (Afield - 1.)*0.50d0

    ! n^a = (d/dv)^a + (A/2) (d/dr)^a
    n_v = 1.
    n_r = (Afield)*0.50d0

    ! s^a = (d/dr)^a
    s_v = 0.
    s_r = 1.

    ! u^a = (d/dv)^a
     u_v = 1.
     u_r = 0.

!get the misner sharp mass compute g^{ab} Sigma,a Sigma,b
       DR = Afield * (1.-xm_reg**2*utau)**2  &
     &      + 2.*(1.-xm_reg**2*utau)*dvusigma

! charge within a sphere, weighted by area element already
! problem as i should be dividing by sqrt(A)!

    output_vec(icharge) = uz * xm_reg**2* &
    &  Sigmafield**2  !* (t_v * Afield - t_r)
 
!http://arxiv.org/abs/1204.4472
    output_vec(imisner) = 0.5*Sigmafield*(Sigmafield**2/rads**2 + 1. - DR)


    ! Charge current components
    J_v = -((q*xm**3.*(-2.*rads**2.*uphir*upii + 2.*rads**2.*uphii*upir + &
         &        2.*q*rads**2.*(uphii**2. + uphir**2.)*uw*xm + &
         &        dxuphir*rads**2.*uphii*xm**3. + &
         &        dxuphir*rads**2.*ualpha*uphii*xm**3. - &
         &        dxuphii*rads**2.*uphir*xm**3. - &
         &        dxuphii*rads**2.*ualpha*uphir*xm**3. + &
         &        dxuphir*uphii*usigma**2.*xm**3. - &
         &        dxuphii*uphir*usigma**2.*xm**3.))/rads**2.)

    J_r = 2.*q*(dxuphir*uphii - dxuphii*uphir)*xm**6.

    ! Dot into n, s vectors
    output_vec(ijn) = n_v*J_v + n_r*J_r
    output_vec(ijs) = s_v*J_v + s_r*J_r

    ! Stress-energy components for phi. (Note T_{\theta\theta},
    ! T_{\phi\phi} also nonzero, but we ignore them here.)
    Tphi_vv = (xm**2.*(4.*rads**4.*upii**2. + 4.*rads**4.*upir**2. -  &
         & 8.*q*rads**4.*uphir*upii*uw*xm + 8.*q*rads**4.*uphii*upir*uw*xm +  &
         & xm**2.*(4.*dxuphii*uphii* &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**3. +  &
         & 4.*dxuphir*uphir*(rads**2. + rads**2.*ualpha + usigma**2.)**2.* &
         & xm**3. + (dxuphii**2. + dxuphir**2.)* &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**4. +  &
         & 4.*uphii**2.*(q**2.*rads**4.*uw**2. +  &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**2.) +  &
         & 4.*uphir**2.*(q**2.*rads**4.*uw**2. +  &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**2.))))/ &
         & (4.*rads**4.)

    Tphi_vr = -((rads**2. + rads**2.*ualpha + usigma**2.)*xm**6.* &
         & (4.*uphii**2. + 4.*uphir**2. + 4.*dxuphii*uphii*xm +  &
         & 4.*dxuphir*uphir*xm + (dxuphii**2. + dxuphir**2.)*xm**2.))/ &
         & (2.*rads**2.) 

    Tphi_rr = xm**6.*(4.*uphii**2. + 4.*uphir**2. + 4.*dxuphii*uphii*xm +  &
         & 4.*dxuphir*uphir*xm + (dxuphii**2. + dxuphir**2.)*xm**2.)

    ! EM
    TA_vv = ((rads**2. + rads**2.*ualpha + usigma**2.)*uz**2.*xm**4.)/(4.*rads**2.)

    TA_vr = -(uz**2.*xm**4.)/4.

    TA_rr = 0.

    ! Dot into vectors
    output_vec(iTphitt) = t_v**2.*Tphi_vv + 2*t_v*t_r*Tphi_vr + t_r**2.*Tphi_rr
    output_vec(iTphits) = t_v*s_v*Tphi_vv + t_v*s_r*Tphi_vr + &
         & t_r*s_v*Tphi_vr + t_r*s_r*Tphi_rr
    output_vec(iTphiss) = s_v**2.*Tphi_vv + 2*s_v*s_r*Tphi_vr + s_r**2.*Tphi_rr

!this is really T_vr - 1/2 g_vr T 
    output_vec(iTphivr) = Tphi_vr


    output_vec(iTAtt) = t_v**2.*TA_vv + 2*t_v*t_r*TA_vr + t_r**2.*TA_rr
    output_vec(iTAts) = t_v*s_v*TA_vv + t_v*s_r*TA_vr + &
         & t_r*s_v*TA_vr + t_r*s_r*TA_rr
    output_vec(iTAss) = s_v**2.*TA_vv + 2*s_v*s_r*TA_vr + s_r**2.*TA_rr

!this is really T_vr - 1/2 g_vr T but T = 0 for the maxwell field
!might want to multiply by 2 to account for R_{AB} g^{AB}
    output_vec(iTAvr) = TA_vr
 
    if (Afield.gt.0.000001) then
       output_vec(iPfluxv) = (0.*TA_vr + Tphi_vr)/sqrt(Afield)
    else
       output_vec(iPfluxv) = 0.0
    end if

    ! Kretschmann
    output_vec(ikretsch) = (xm**4.*(4.*usigma**4. + 16.*dxus**2.*rads**4.*xm**2. -  &
         & 16.*rads**4.*upii**2.*utau*xm**2. - 16.*rads**4.*upir**2.*utau*xm**2. +  &
         & 32.*dxus*rads**2.*usigma*utau*xm**2. -  &
         & 8.*rads**2.*usigma**2.*utau*xm**2. - 8.*usigma**4.*utau*xm**2. +  &
         & 32.*usigma**2.*utau**2.*xm**2. - 8.*dxutau*rads**4.*upii**2.*xm**3. -  &
         & 8.*dxutau*rads**4.*upir**2.*xm**3. +  &
         & 32.*dxus**2.*rads**4.*usigma*xm**3. +  &
         & 8.*dxubeta*rads**2.*usigma*utau*xm**3. -  &
         & 48.*rads**4.*upii**2.*usigma*utau*xm**3. -  &
         & 48.*rads**4.*upir**2.*usigma*utau*xm**3. +  &
         & 16.*dxutau*usigma**2.*utau*xm**3. +  &
         & 64.*dxus*rads**2.*usigma**2.*utau*xm**3. +  &
         & 96.*usigma**3.*utau**2.*xm**3. + 16.*usigma*utau**3.*xm**3. +  &
         & 32.*q*rads**4.*uphir*upii*utau*uw*xm**3. -  &
         & 32.*q*rads**4.*uphii*upir*utau*uw*xm**3. +  &
         & dxubeta**2.*rads**4.*xm**4. +  &
         & 4.*dxubeta*dxutau*rads**2.*usigma*xm**4. -  &
         & 24.*dxutau*rads**4.*upii**2.*usigma*xm**4. -  &
         & 24.*dxutau*rads**4.*upir**2.*usigma*xm**4. +  &
         & 4.*dxutau**2.*usigma**2.*xm**4. +  &
         & 16.*dxus**2.*rads**4.*usigma**2.*xm**4. +  &
         & 32.*dxubeta*rads**2.*usigma**2.*utau*xm**4. -  &
         & 48.*rads**4.*upii**2.*usigma**2.*utau*xm**4. -  &
         & 48.*rads**4.*upir**2.*usigma**2.*utau*xm**4. +  &
         & 64.*dxutau*usigma**3.*utau*xm**4. +  &
         & 32.*dxus*rads**2.*usigma**3.*utau*xm**4. +  &
         & 4.*dxubeta*rads**2.*utau**2.*xm**4. + 4.*rads**4.*utau**2.*xm**4. +  &
         & 8.*dxutau*usigma*utau**2.*xm**4. +  &
         & 8.*rads**2.*usigma**2.*utau**2.*xm**4. +  &
         & 116.*usigma**4.*utau**2.*xm**4. + 64.*usigma**2.*utau**3.*xm**4. +  &
         & 4.*utau**4.*xm**4. + 16.*dxutau*q*rads**4.*uphir*upii*uw*xm**4. -  &
         & 16.*dxutau*q*rads**4.*uphii*upir*uw*xm**4. +  &
         & 96.*q*rads**4.*uphir*upii*usigma*utau*uw*xm**4. -  &
         & 96.*q*rads**4.*uphii*upir*usigma*utau*uw*xm**4. -  &
         & 16.*q**2.*rads**4.*uphii**2.*utau*uw**2.*xm**4. -  &
         & 16.*q**2.*rads**4.*uphir**2.*utau*uw**2.*xm**4. +  &
         & 4.*dxubeta**2.*rads**4.*usigma*xm**5. +  &
         & 16.*dxubeta*dxutau*rads**2.*usigma**2.*xm**5. -  &
         & 24.*dxutau*rads**4.*upii**2.*usigma**2.*xm**5. -  &
         & 24.*dxutau*rads**4.*upir**2.*usigma**2.*xm**5. +  &
         & 16.*dxutau**2.*usigma**3.*xm**5. +  &
         & 48.*dxubeta*rads**2.*usigma**3.*utau*xm**5. -  &
         & 16.*rads**4.*upii**2.*usigma**3.*utau*xm**5. -  &
         & 16.*rads**4.*upir**2.*usigma**3.*utau*xm**5. +  &
         & 96.*dxutau*usigma**4.*utau*xm**5. +  &
         & 16.*dxubeta*rads**2.*usigma*utau**2.*xm**5. +  &
         & 32.*dxutau*usigma**2.*utau**2.*xm**5. +  &
         & 64.*usigma**5.*utau**2.*xm**5. + 96.*usigma**3.*utau**3.*xm**5. +  &
         & 16.*usigma*utau**4.*xm**5. +  &
         & 48.*dxutau*q*rads**4.*uphir*upii*usigma*uw*xm**5. -  &
         & 48.*dxutau*q*rads**4.*uphii*upir*usigma*uw*xm**5. +  &
         & 96.*q*rads**4.*uphir*upii*usigma**2.*utau*uw*xm**5. -  &
         & 96.*q*rads**4.*uphii*upir*usigma**2.*utau*uw*xm**5. -  &
         & 8.*dxutau*q**2.*rads**4.*uphii**2.*uw**2.*xm**5. -  &
         & 8.*dxutau*q**2.*rads**4.*uphir**2.*uw**2.*xm**5. -  &
         & 48.*q**2.*rads**4.*uphii**2.*usigma*utau*uw**2.*xm**5. -  &
         & 48.*q**2.*rads**4.*uphir**2.*usigma*utau*uw**2.*xm**5. +  &
         & 6.*dxubeta**2.*rads**4.*usigma**2.*xm**6. +  &
         & 24.*dxubeta*dxutau*rads**2.*usigma**3.*xm**6. -  &
         & 8.*dxutau*rads**4.*upii**2.*usigma**3.*xm**6. -  &
         & 8.*dxutau*rads**4.*upir**2.*usigma**3.*xm**6. +  &
         & 24.*dxutau**2.*usigma**4.*xm**6. +  &
         & 32.*dxubeta*rads**2.*usigma**4.*utau*xm**6. +  &
         & 64.*dxutau*usigma**5.*utau*xm**6. +  &
         & 24.*dxubeta*rads**2.*usigma**2.*utau**2.*xm**6. +  &
         & 48.*dxutau*usigma**3.*utau**2.*xm**6. +  &
         & 16.*usigma**6.*utau**2.*xm**6. + 64.*usigma**4.*utau**3.*xm**6. +  &
         & 24.*usigma**2.*utau**4.*xm**6. +  &
         & 48.*dxutau*q*rads**4.*uphir*upii*usigma**2.*uw*xm**6. -  &
         & 48.*dxutau*q*rads**4.*uphii*upir*usigma**2.*uw*xm**6. +  &
         & 32.*q*rads**4.*uphir*upii*usigma**3.*utau*uw*xm**6. -  &
         & 32.*q*rads**4.*uphii*upir*usigma**3.*utau*uw*xm**6. -  &
         & 24.*dxutau*q**2.*rads**4.*uphii**2.*usigma*uw**2.*xm**6. -  &
         & 24.*dxutau*q**2.*rads**4.*uphir**2.*usigma*uw**2.*xm**6. -  &
         & 48.*q**2.*rads**4.*uphii**2.*usigma**2.*utau*uw**2.*xm**6. -  &
         & 48.*q**2.*rads**4.*uphir**2.*usigma**2.*utau*uw**2.*xm**6. +  &
         & 4.*dxubeta**2.*rads**4.*usigma**3.*xm**7. +  &
         & 16.*dxubeta*dxutau*rads**2.*usigma**4.*xm**7. +  &
         & 16.*dxutau**2.*usigma**5.*xm**7. +  &
         & 8.*dxubeta*rads**2.*usigma**5.*utau*xm**7. +  &
         & 16.*dxutau*usigma**6.*utau*xm**7. +  &
         & 16.*dxubeta*rads**2.*usigma**3.*utau**2.*xm**7. +  &
         & 32.*dxutau*usigma**4.*utau**2.*xm**7. +  &
         & 16.*usigma**5.*utau**3.*xm**7. + 16.*usigma**3.*utau**4.*xm**7. +  &
         & 16.*dxutau*q*rads**4.*uphir*upii*usigma**3.*uw*xm**7. -  &
         & 16.*dxutau*q*rads**4.*uphii*upir*usigma**3.*uw*xm**7. -  &
         & 24.*dxutau*q**2.*rads**4.*uphii**2.*usigma**2.*uw**2.*xm**7. -  &
         & 24.*dxutau*q**2.*rads**4.*uphir**2.*usigma**2.*uw**2.*xm**7. -  &
         & 16.*q**2.*rads**4.*uphii**2.*usigma**3.*utau*uw**2.*xm**7. -  &
         & 16.*q**2.*rads**4.*uphir**2.*usigma**3.*utau*uw**2.*xm**7. +  &
         & dxubeta**2.*rads**4.*usigma**4.*xm**8. +  &
         & 4.*dxubeta*dxutau*rads**2.*usigma**5.*xm**8. +  &
         & 4.*dxutau**2.*usigma**6.*xm**8. +  &
         & 4.*dxubeta*rads**2.*usigma**4.*utau**2.*xm**8. +  &
         & 8.*dxutau*usigma**5.*utau**2.*xm**8. + 4.*usigma**4.*utau**4.*xm**8. -  &
         & 8.*dxutau*q**2.*rads**4.*uphii**2.*usigma**3.*uw**2.*xm**8. -  &
         & 8.*dxutau*q**2.*rads**4.*uphir**2.*usigma**3.*uw**2.*xm**8. +  &
         & 4.*ubeta**2.*xm**2.*(rads + rads*usigma*xm)**4. +  &
         & 4.*rads**2.*ubeta*xm**2.*(1. + usigma*xm)**4.* &
         & (4.*usigma*utau + dxubeta*rads**2.*xm + 2.*dxutau*usigma*xm +  &
         & 2.*utau**2.*xm) + 16.*rads**4.*us**2.*(-1. + utau*xm**2.)**2. +  &
         & 16.*rads**2.*us*(-1. + utau*xm**2.)* &
         & (rads**2.*utau*xm**2. + usigma**2.*(-1. + utau*xm**2.))))/ &
         & (rads + rads*usigma*xm)**4.

    ! \rho^2 d_+\Sigma
    output_vec(idplussigma) = us*xm**2. + 0.5*(xm*usigma + 1.)**2./rads + 0.5*xm**2.

    ! \rho * d_+d_+\Sigma equation (an independent residual)
    output_vec(ires1) = (4.*dvusigma - 2.*dxus*xm + 4.*dvus*rads**2.*xm + ubeta*xm + &
         & 4.*dvusigma*usigma*xm + 2.*rads**2.*upii**2.*xm**2. + &
         & 2.*rads**2.*upir**2.*xm**2. - 4.*dxus*usigma*xm**2. + &
         & 2.*ubeta*usigma*xm**2. - 2.*dxus*rads**2.*xm**3. + &
         & rads**2.*ubeta*xm**3. + 2.*rads**2.*upii**2.*usigma*xm**3. + &
         & 2.*rads**2.*upir**2.*usigma*xm**3. - 2.*dxus*usigma**2.*xm**3. + &
         & ubeta*usigma**2.*xm**3. - 4.*q*rads**2.*uphir*upii*uw*xm**3. + &
         & 4.*q*rads**2.*uphii*upir*uw*xm**3. - &
         & 4.*q*rads**2.*uphir*upii*usigma*uw*xm**4. + &
         & 4.*q*rads**2.*uphii*upir*usigma*uw*xm**4. + &
         & 2.*q**2.*rads**2.*uphii**2.*uw**2.*xm**4. + &
         & 2.*q**2.*rads**2.*uphir**2.*uw**2.*xm**4. + &
         & 2.*q**2.*rads**2.*uphii**2.*usigma*uw**2.*xm**5. + &
         & 2.*q**2.*rads**2.*uphir**2.*usigma*uw**2.*xm**5. - &
         & 2.*ualpha*(-1. + utau*xm**2. + dxus*rads**2.*xm**3. + &
         & usigma*xm*(-1. + utau*xm**2.)) + &
         & 2.*us*(rads**2.*ubeta*xm**3. + &
         & 2.*(1. + usigma*xm)*(-1. + utau*xm**2.)))/(4.*rads**2.)

    ! (d_+W)' equation (another independent residual)
    output_vec(ires2) = (xm**2.*(-4.*rads**2.*us*uz*xm + &
         & (1. + usigma*xm)*(dxuz - 2.*dvuz*rads**2. + &
         & 4.*q*rads**2.*uphir*upii*xm - &
         & 4.*q*rads**2.*uphii*upir*xm + 2.*dxuz*usigma*xm + &
         & dxuz*rads**2.*xm**2. + dxuz*rads**2.*ualpha*xm**2. + &
         & dxuz*usigma**2.*xm**2. - &
         & 4.*q**2.*rads**2.*(uphii**2. + uphir**2.)*uw*xm**2.) + &
         & 2.*uz*(rads**2.*ualpha*xm*(1. + usigma*xm) + &
         & usigma*(1. + 2.*usigma*xm + rads**2.*xm**2. + &
         & usigma**2.*xm**2.))))/(2.*rads**2.*(1. + usigma*xm))


  end subroutine aux_output_rhs

  subroutine horizon_output_rhs(u,dxu,xm,dx,out_vec)
    use m_pars, only : isigma, itau, is, ialpha, ibeta, iw, iz, ipir, &
         & ipii, iphii, iphir, ihorjs, ihorTphits, ihorTAts, q, rads,ihorcharge, N
    implicit none
    real*8, dimension(:), intent(in) :: u, dxu
    real*8, intent(in) :: xm, dx
    real*8, dimension(:), intent(inout) :: out_vec
    real*8 :: uphir,uphii,usigma,utau,us,ualpha,ubeta,uw,uz,upir,upii
    real*8 :: dxuphir,dxuphii,dxutau,dxus,dxubeta,dxuz
    real*8 :: Tphi_vr, Tphi_vv, Tphi_rr, TA_vr, TA_vv, TA_rr, J_v, J_r
    real*8 :: t_v, t_r, n_v, n_r, s_v, s_r, xm_reg

    real*8 :: Afield, Sigmafield



    ! Needed variables
    uphir = u(iphir)
    dxuphir = dxu(iphir)

    uphii = u(iphii)
    dxuphii = dxu(iphii)    

    usigma = u(isigma)

    utau = u(itau)
    dxutau = dxu(itau)

    us = u(is)
    dxus = dxu(is)

    ualpha = u(ialpha)

    ubeta = u(ibeta)
    dxubeta = dxu(ibeta)

    uw = u(iw)

    uz = u(iz)
    dxuz = dxu(iz)

    upir = u(ipir)
    upii = u(ipii)

    ! Define vectors (up-index):

    ! Use a regulator, since A blows up at infinity (ad hoc). If
    ! fields go smoothly to 0 in the end, this is fine.
    if (xm.lt.0.0001*dx) then
       xm_reg = 0.0001*dx
    else
       xm_reg = xm
    end if

    ! t^a = (d/dv)^a + (A-1)/2 (d/dr)^a
    t_v = 1.
    t_r = (ualpha + (usigma + 1./xm_reg)**2./rads**2.) *0.5d0

    ! n^a = (d/dv)^a + (A/2) (d/dr)^a
    n_v = 1.
    n_r = (ualpha + (usigma + 1./xm_reg)**2./rads**2. + 1.) *0.5d0

    ! s^a = (d/dr)^a
    s_v = 0.
    s_r = 1.


! charge within a sphere, weighted by area element already
! problem as i should be dividing by sqrt(A)!
! this is compensated by area element?

    Afield = (ualpha + (usigma + 1./xm_reg)**2./rads**2. + 1.)
    Sigmafield = usigma + 1./xm

    !out_vec(ihorcharge) = - uz * xm**2* &
    !&  Sigmafield**2  * (t_v * Afield - t_r)

    out_vec(ihorcharge) =  uz * xm**2* Sigmafield**2


    ! Charge current components
    J_v = -((q*xm**3.*(-2.*rads**2.*uphir*upii + 2.*rads**2.*uphii*upir + &
         &        2.*q*rads**2.*(uphii**2. + uphir**2.)*uw*xm + &
         &        dxuphir*rads**2.*uphii*xm**3. + &
         &        dxuphir*rads**2.*ualpha*uphii*xm**3. - &
         &        dxuphii*rads**2.*uphir*xm**3. - &
         &        dxuphii*rads**2.*ualpha*uphir*xm**3. + &
         &        dxuphir*uphii*usigma**2.*xm**3. - &
         &        dxuphii*uphir*usigma**2.*xm**3.))/rads**2.)

    J_r = 2.*q*(dxuphir*uphii - dxuphii*uphir)*xm**6.

    ! Dot into s vector
    !out_vec(ijn) = n_v*J_v + n_r*J_r
    out_vec(ihorjs) = s_v*J_v + s_r*J_r

    ! Stress-energy components for phi. (Note T_{\theta\theta},
    ! T_{\phi\phi} also nonzero, but we ignore them here.)
    Tphi_vv = (xm**2.*(4.*rads**4.*upii**2. + 4.*rads**4.*upir**2. -  &
         & 8.*q*rads**4.*uphir*upii*uw*xm + 8.*q*rads**4.*uphii*upir*uw*xm +  &
         & xm**2.*(4.*dxuphii*uphii* &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**3. +  &
         & 4.*dxuphir*uphir*(rads**2. + rads**2.*ualpha + usigma**2.)**2.* &
         & xm**3. + (dxuphii**2. + dxuphir**2.)* &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**4. +  &
         & 4.*uphii**2.*(q**2.*rads**4.*uw**2. +  &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**2.) +  &
         & 4.*uphir**2.*(q**2.*rads**4.*uw**2. +  &
         & (rads**2. + rads**2.*ualpha + usigma**2.)**2.*xm**2.))))/ &
         & (4.*rads**4.)

    Tphi_vr = -((rads**2. + rads**2.*ualpha + usigma**2.)*xm**6.* &
         & (4.*uphii**2. + 4.*uphir**2. + 4.*dxuphii*uphii*xm +  &
         & 4.*dxuphir*uphir*xm + (dxuphii**2. + dxuphir**2.)*xm**2.))/ &
         & (2.*rads**2.) 

    Tphi_rr = xm**6.*(4.*uphii**2. + 4.*uphir**2. + 4.*dxuphii*uphii*xm +  &
         & 4.*dxuphir*uphir*xm + (dxuphii**2. + dxuphir**2.)*xm**2.)

    ! EM
    TA_vv = ((rads**2. + rads**2.*ualpha + usigma**2.)*uz**2.*xm**4.)/(4.*rads**2.)

    TA_vr = -(uz**2.*xm**4.)/4.

    TA_rr = 0.

    ! Dot into vectors
    !out_vec(ihorTphitt) = t_v**2.*Tphi_vv + 2*t_v*t_r*Tphi_vr + t_r**2.*Tphi_rr
    out_vec(ihorTphits) = t_v*s_v*Tphi_vv + t_v*s_r*Tphi_vr + &
         & t_r*s_v*Tphi_vr + t_r*s_r*Tphi_rr

    !out_vec(ihorTAtt) = t_v**2.*TA_vv + 2*t_v*t_r*TA_vr + t_r**2.*TA_rr
    out_vec(ihorTAts) = t_v*s_v*TA_vv + t_v*s_r*TA_vr + &
         & t_r*s_v*TA_vr + t_r*s_r*TA_rr
 

  end subroutine horizon_output_rhs

end module M_RHS
