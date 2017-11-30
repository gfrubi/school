!this file contains the RHSs of the equations.
!divided in those needed for the time integration
!and those for the space integration

MODULE M_RHS
  implicit none

contains


  subroutine rhs_time(u,du,source,dx,N,x,time)
    use m_pars, only : amp,IMAX, isigma, itau, is, ialpha, ibeta, ipi, iphi, &
         & lambda, p2, a2, madm, rhoads, rads, scalaronly
    implicit none
    real*8, dimension(:,:), intent(in) :: u, du
    real*8, dimension(:,:), intent(inout) :: source
    real*8, dimension(:), intent(in) :: x
    real*8 :: dx, v,time, uphi,usigma,utau,us,ualpha,ubeta,upi,dxuphi
    integer :: i, N
    logical :: ltrace
    real*8 :: tangh,p0,dt_p0,d2t_p0,d3t_p0, xm, Axx, aa

    !for debug
    ltrace=.false.

    !here in case we want to keep messy things in 
    !a separate module	
    source(:,:) = 0.0

    ! Boundary evolution
    source(iphi,1) = 0.0d0

    ! Evolution away from boundary
    do i=2, IMAX
       if (scalaronly) then
          ! Metric function A
          aa = 1 - 2.*madm*x(i) + (rhoads**2.) / (x(i)**2.)
          source(iphi,i) = u(ipi,i)/x(i) + aa*x(i)*u(iphi,i) &
               & + 0.5*aa*(x(i)**2.)*du(iphi,i)
       else
          uphi = u(iphi,i)
          dxuphi = du(iphi,i)
          usigma = u(isigma,i)
          utau = u(itau,i)
          us = u(is,i)
          ualpha = u(ialpha,i)
          ubeta = u(ibeta,i)
          upi = u(ipi,i)
          
          source(iphi,i) = (2.*rads**2.*upi + (2.*uphi + dxuphi*x(i))* &
               & (1. + 2.*usigma*x(i) + rads**2.*x(i)**2. + &
               & rads**2.*ualpha*x(i)**2. + usigma**2.*x(i)**2.))/ &
               & (2.*rads**2.*x(i))
       end if
    end do
    !     source(iphi,1) = 0d0 !.0.5d0*u(ig,1)
    !	source(ig,IMAX) = -du(ig,IMAX)
    !	print*,'source',source(iphi,1)
    !	source(ig,IMAX) = -0.5d0*(du(iphi,IMAX)-du(iphi,IMAX-1))/dx


  end subroutine rhs_time

  subroutine rhs_space(uguess,u,du,d2u,source,dx,Nv,Nspatial,N,xm,time)
    use m_pars, only : amp, gamma,isigma, itau, is, ialpha, ibeta, ipi, iphi,lambda,a2,p2,rhoads,madm,rads,scalaronly
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
    real*8 :: uphi,usigma,utau,us,ualpha,ubeta,upi,bet,bb,alp,aa
    real*8 :: dx_ppi,dxuphi,dx_bet,dx_bb,dx_alp,dx_aa
    real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,t0,rho
    real*8 :: RR(3), lnxm
    !for debug
    ltrace=.false.


    !define first all that is needed, notice phi is evolved
    !already so they come from u, the rest from uguess

    uphi = u(iphi)
    dxuphi = du(iphi)

    usigma = uguess(isigma)
    utau = uguess(itau)
    us = uguess(is)
    ualpha = uguess(ialpha)
    ubeta = uguess(ibeta)
    upi = uguess(ipi)

    source = 0.

    ! First calculate source for scalar field
    
    ! We have a special case at rho = 0 end point
    if (xm.lt.(0.1*dx)) then
       source(ipi) = - 1.5d0*(rhoads**2.)*dxuphi
    else
       if (scalaronly) then
          ! Metric function A corresponding to Schwarzschild-AdS
          aa = 1.d0 - 2.0d0*madm*xm + (rhoads**2.) / (xm**2.)
          source(ipi) = - aa*xm*uphi - 0.5d0*aa*(xm**2.)*dxuphi
       else
          source(ipi) = -(2.*rads**2.*upi*xm*(usigma + utau*xm) +  &
               &  (2.*uphi + dxuphi*xm)* &
               &  (1. + 2.*usigma*xm + rads**2.*xm**2. + 2.*rads**2.*us*xm**2. +  &
               &  usigma**2.*xm**2.))/(2.*rads**2.*xm*(1. + usigma*xm))
       end if
    end if

    if (.not.scalaronly) then
       ! End point
       if (xm.lt.(0.1*dx)) then
          source(is) = -madm
          source(ialpha) = -2.*madm
          source(ibeta) = 0.
          source(isigma) = 0.
          source(itau) = 0.
       else
          source(is) = -(utau*(3. + 6.*usigma*xm + rads**2.*xm**2. +  &
               &  3.*usigma**2.*xm**2.) + 2.*rads**2.*us*(-1. + utau*xm**2.))/ &
               &  (2.*rads**2.*xm*(1. + usigma*xm))
          
          source(ialpha) = ubeta
          
          source(ibeta) = (2.*utau + 4.*usigma*utau*xm + 4.*uphi**2.*xm**2. +  &
               &  4.*rads**2.*uphi*upi*xm**2. - 2.*rads**2.*utau*xm**2. +  &
               &  2.*usigma**2.*utau*xm**2. - 2.*utau**2.*xm**2. +  &
               &  4.*dxuphi*uphi*xm**3. + 2.*dxuphi*rads**2.*upi*xm**3. +  &
               &  16.*uphi**2.*usigma*xm**3. +  &
               &  8.*rads**2.*uphi*upi*usigma*xm**3. -  &
               &  4.*usigma*utau**2.*xm**3. + dxuphi**2.*xm**4. +  &
               &  16.*dxuphi*uphi*usigma*xm**4. +  &
               &  4.*dxuphi*rads**2.*upi*usigma*xm**4. +  &
               &  24.*uphi**2.*usigma**2.*xm**4. +  &
               &  4.*rads**2.*uphi*upi*usigma**2.*xm**4. -  &
               &  2.*usigma**2.*utau**2.*xm**4. + 4.*dxuphi**2.*usigma*xm**5. +  &
               &  24.*dxuphi*uphi*usigma**2.*xm**5. +  &
               &  2.*dxuphi*rads**2.*upi*usigma**2.*xm**5. +  &
               &  16.*uphi**2.*usigma**3.*xm**5. +  &
               &  6.*dxuphi**2.*usigma**2.*xm**6. +  &
               &  16.*dxuphi*uphi*usigma**3.*xm**6. +  &
               &  4.*uphi**2.*usigma**4.*xm**6. + 4.*dxuphi**2.*usigma**3.*xm**7. +  &
               &  4.*dxuphi*uphi*usigma**4.*xm**7. +  &
               &  dxuphi**2.*usigma**4.*xm**8. -  &
               &  2.*ubeta*xm*(rads + rads*usigma*xm)**2. -  &
               &  4.*rads**2.*us*(-1. + utau*xm**2.))/ &
               &  (rads**2.*xm**2.*(1. + usigma*xm)**2.)
     
          source(isigma) = utau

          source(itau) = -(4.*utau + xm**2.*(2.*uphi + dxuphi*xm)**2.*(1. + usigma*xm))/(2.*xm)

       end if
    end if

  end subroutine rhs_space


end module M_RHS
