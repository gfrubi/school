!this file contains the RHSs of the equations.
!divided in those needed for the time integration
!and those for the space integration

MODULE M_RHS
  implicit none

contains


  subroutine rhs_time(u,du,source,dx,N,x,time)
    use m_pars, only : amp,IMAX,iphi, ipi, ib,ibeta,ia, ialp, &
         & lambda, p2, a2, mbh, rhoads
    implicit none
    real*8, dimension(:,:), intent(in) :: u, du
    real*8, dimension(:,:), intent(inout) :: source
    real*8, dimension(:), intent(in) :: x
    real*8 :: dx, v,time
    integer :: i, N
    logical :: ltrace
    real*8 :: tangh,p0,dt_p0,d2t_p0,d3t_p0, xm, Axx, aa

    !for debug
    ltrace=.false.

    !here in case we want to keep messy things in 
    !a separate module	
    source(:,:) = 0.0

    ! Boundary evolution
    source(iphi,1) = 0.0d0 ! du(ipi,1) + 1.5d0*(rhoads**2.)*du(iphi,1)

    ! Evolution away from boundary
    do i=2, IMAX
       ! Metric function A
       aa = 1 - 2.*mbh*x(i) + (rhoads**2.) / (x(i)**2.)
       source(iphi,i) = u(ipi,i)/x(i) + aa*x(i)*u(iphi,i) &
            & + 0.5*aa*(x(i)**2.)*du(iphi,i)
    end do
    !     source(iphi,1) = 0d0 !.0.5d0*u(ig,1)
    !	source(ig,IMAX) = -du(ig,IMAX)
    !	print*,'source',source(iphi,1)
    !	source(ig,IMAX) = -0.5d0*(du(iphi,IMAX)-du(iphi,IMAX-1))/dx


  end subroutine rhs_time

  subroutine rhs_space(uguess,u,du,d2u,source,dx,N,xm,time)
    use m_pars, only : amp, gamma,iphi, ipi,ibeta,ia, ialp,lambda,a2,p2,rhoads,mbh
    implicit none
    real*8, dimension(:), intent(in) :: u, du
    real*8 :: time, xm

    real*8, dimension(6) :: source, d2u
    real*8, dimension(6) :: uguess
    real*8 :: dx, v
    integer :: i, j, N
    logical :: ltrace

    real*8, dimension(6) :: At,Arr,drArr,Ar,A0,J0
    real*8 :: tangh,p0,dt_p0,d2t_p0,d3t_p0
    real*8 :: uppi,pph,bet,bb,alp,aa
    real*8 :: dx_ppi,dx_pph,dx_bet,dx_bb,dx_alp,dx_aa
    real*8 :: s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,t0,rho
    real*8 :: RR(3), lnxm
    !for debug
    ltrace=.false.


    !define first all that is needed, notice phi and b are evolved
    !already so they come from u, the rest from uguess

    pph = u(iphi)
    uppi = uguess(ipi)

    dx_pph = du(iphi)

    ! We have a special case at rho = 0 end point
    if (xm.lt.(0.25*dx)) then
       source(ipi) = - 1.5d0*(rhoads**2.)*dx_pph
    else
       ! Metric function A
       aa = 1.d0 - 2.0d0*mbh*xm + (rhoads**2.) / (xm**2.)
       source(ipi) = - aa*xm*pph - 0.5d0*aa*(xm**2.)*dx_pph
    end if

    !print *, source(ipi)

  end subroutine rhs_space


end module M_RHS
