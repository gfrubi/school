MODULE M_INITIAL
  use M_pars
  use M_derivs
  implicit none

CONTAINS

  subroutine initial(Uold,dx,x,Nv,N,dt,lambda)
    use m_pars, only : iphir, iphii, disip, amp, rhoin, rhoout
    implicit none
    real*8, dimension(:,:) :: Uold
    real*8, dimension(:) :: x
    real*8 :: dx, dy,dt,lambda,xmin,xmax
    integer :: Nv, N, i, ret, gft_out_brief

    real*8 :: ppi
    real*8 :: p0,Asx,As,xbc,lnxbc,s,sl,sr

    ppi = 3.141592653589793
    p0 = amp

    !initialize all to zero, then fill
    uold(:,:) = 0.0d0
    !

    !giving data for the field g, the field phi
    !will have to be integrated radially to be obtained
    !this we do as part of the main routine

    xmin=rhoin + (rhoout-rhoin)*0.25
    xmax=rhoin + (rhoout-rhoin)*0.55
    
    do i=1,N
       if(x(i).ge.xmin.and.x(i).lt.xmax) then
          uold(iphir,i) = -amp *(x(i)-xmin)**3*(x(i)-xmax)**3*(1.+sin(x(i)*100./(rhoin-rhoout)))
       end if

!!$       if(x(i).ge.1.85.and.x(i).lt.1.95) then
!!$          uold(ig,i) = 1e10*amp *(x(i)-1.85)**3*(x(i)-1.95)**3*(2.*x(i)-1.85-1.95)
!!$       end if
!!$
    end do

    print *, "Set up initial data."

  end subroutine initial


END MODULE M_INITIAL
