!this file contains the RHSs of the equations.
!divided in those needed for the time integration
!and those for the space integration

MODULE M_RHS
  implicit none

contains


  subroutine rhs_time(u,du,source,dx,N,x,time)
    use m_pars, only : amp,IMAX,iphi, ig, ib,ibeta,ia, ialp, lambda, p2, a2
    implicit none
    real*8, dimension(:,:), intent(in) :: u, du
    real*8, dimension(:,:), intent(inout) :: source
    real*8, dimension(:), intent(in) :: x
    real*8 :: dx, v,time
    integer :: i, N
    logical :: ltrace

    !for debug
    ltrace=.false.

    !here in case we want to keep messy things in 
    !a separate module	
    source(:,:) = 0.0

    do i=1, IMAX-1
       source(ig,i) = du(ig,i) * 0.5d0 
    end do
    source(iphi,1) = 0d0 

  end subroutine rhs_time


  subroutine rhs_space(uguess,u,du,d2u,source,dx,N,xm,time)
    use m_pars, only : amp, gamma,iphi, ig,ibeta,ia, ialp,lambda,a2,p2
    implicit none
    real*8, dimension(:), intent(in) :: u, du
    real*8 :: time, xm

    real*8, dimension(6) :: source, d2u
    real*8, dimension(6) :: uguess
    real*8 :: dx, v
    integer :: N
    logical :: ltrace

    real*8 :: gg
    real*8 :: pph
    
    !for debug
    ltrace=.false.

    pph = u(iphi)
    gg = u(ig)

    source(iphi)  = gg

  end subroutine rhs_space


end module M_RHS
