!this gives and error that goes ~ dt^2 , dx^4, so if 
!the spacing is sufficiently small, one sees the dx^4 error...

module M_EVOLVE
  implicit none

contains

  ! try getting higher order, rk in time and space in 
  ! sequential order 

  subroutine evolve(unew,uold,dx,dt,time,x,Nv,N,lambda)
    use m_derivs
    use m_rhs
    use m_pars, only :iphi, ig, ib,ibeta,ia, ialp, disip, &
         &     a2iold, a2i, p2, p2old, amp, IMAX, a2, mask, freqhor,initime,&
         &     area_rad, IMAXOLD, null_rad, ITE

    implicit none

    real*8, dimension(:,:), intent(inout):: unew
    real*8, dimension(:,:), intent(in):: uold
    real*8, dimension(:):: x	
    real*8 dx, dt, time,lambda
    integer :: Nv,N

    !here are some basic local variables, modify as needed

    real*8 :: xm
    real*8, dimension(:,:), allocatable :: tempdis, u1,u2,u3,u4,du
    real*8, dimension(:,:), allocatable :: urk1,urk2,urk3,urk4
    real*8 :: ak1, ak2, ak3, ak4, akt

    integer i,j
    logical :: ltrace

    !for debug
    ltrace = .false.
    !allocate memory
    allocate(tempdis(Nv,N),u1(Nv,N),u2(Nv,N),u3(Nv,N),u4(Nv,N),du(Nv,N), &
         &           urk1(Nv,N),urk2(Nv,N),urk3(Nv,N),urk4(Nv,N))

    !integrate here so that the fields are field with proper values initially
    if(time.lt.initime+0.25*dt) then
       call boundary_fix(uold,x,dx,time,dt,ak1,0)
       call rk_fullradial(uold,du,x,dx,dt,N,0,time)	    
    end if
    !LUIS
    u1 = uold
    u2 = uold
    u3 = uold
    u4 = uold

    !march up 1st step rk, then full radial out
    !notice the rk evoln step is always the same, only changing
    !by the size of the step, which we pass to that routine
    call rk_timeup(uold,uold,u1,du,tempdis,urk1,x,dx,0.5*dt,N,time)
    call boundary_fix(u1,x,dx,time,dt*0.5,ak1,0)
    call rk_fullradial(u1,du,x,dx,dt,N,0,time+0.5*dt)

    !now do the 2nd RK step, then full radial out
    call rk_timeup(uold,u1,u2,du,tempdis,urk2,x,dx,0.5*dt,N,time+0.5*dt)
    call boundary_fix(u2,x,dx,time+0.5*dt,dt*0.5,ak2,0)
    call rk_fullradial(u2,du,x,dx,dt,N,0,time+0.5*dt)


    !now do the 3rd RK step, then full radial out
    call rk_timeup(uold,u2,u3,du,tempdis,urk3,x,dx,dt,N,time+dt)
    call boundary_fix(u3,x,dx,time+dt,dt,ak3,0)
    call rk_fullradial(u3,du,x,dx,dt,N,0,time+dt)


    !now do the 4th RK step, then full radial ou
    call rk_timeup(uold,u3,u4,du,tempdis,urk4,x,dx,dt,N,time+dt)
    call boundary_fix(u4,x,dx,time+dt,dt,ak4,0)

    !now combine the RKs
    unew(ig,1:IMAX)=uold(ig,1:IMAX)+dt/6.*(urk1(ig,1:IMAX)+2.*urk2(ig,1:IMAX)+ &
         &                                       2.*urk3(ig,1:IMAX)+urk4(ig,1:IMAX))
    unew(ig,IMAX)=unew(ig,1)


    call rk_fullradial(unew,du,x,dx,dt,N,0,time+dt)


    deallocate(tempdis,u1,u2,u3,u4,du,urk1,urk2,urk3,urk4)

    return

  end subroutine evolve

  ! RK intermediate step for the time advance
  subroutine rk_timeup(uold,uprev,unew,du,tempdis,urk,x,dx,dt,N,time)
    use m_derivs
    use m_rhs
    use m_pars, only : iphi, ipi, ig,ibeta,ia, ialp, disip, IMAX, mask
    implicit none
    real*8, dimension(:,:) :: uold, uprev,unew, du,urk,tempdis
    real*8, dimension(:) :: x
    real*8 :: dx, dt, time
    integer :: N

    integer :: i, j

    do j=1,2
       call derivs(uprev(j,:),du(j,:),dx, N,1,IMAX)
       call dissip(uprev(j,:),tempdis(j,:),dx, N,1,IMAX)
    end do

    call rhs_time(uprev,du,urk,dx,N,x,time)

    do i = 1, IMAX         
       unew(ig,i) = uold(ig,i)+ dt*(urk(ig,i)+disip*tempdis(ig,i))
    end do

  end subroutine rk_timeup

  !RK full march radially
  subroutine rk_fullradial(unew,du,x,dx,dt,N,pp,time)
    use m_derivs
    use m_rhs
    use m_pars, only :iphi, ipi, ib,ibeta,ia, ialp, disip, &
         &                    a2iold, a2i, p2, p2old, amp,lambda, IMAX, a2
    implicit none
    real*8, dimension(:,:) :: unew, du
    real*8, dimension(:) :: x
    real*8 :: dx, dt,time
    integer :: N, pp, ll 

    real*8 :: xm,del
    real*8, dimension(6) :: k1,k2,k3,k4, u1,u2,u3,u4, dum,uguess, um
    integer :: i,j,I0

    integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
    real*8, dimension(:,:), allocatable :: d2u
    real*8 :: coeff(3), p0, dtp0, d2tp0, pixx

    allocate(d2u(6,N))


    !below is a call to derivatives that would be needed
    !if the rhs of the radial equation has them 
    !in this toy model they are not needed, 
    !but will leave here in case they are needed
    !in the future      
    do j=1,6
       call derivs(unew(j,:),du(j,:),dx, N,1,IMAX)
       call dersecond2(unew(j,:),d2u(j,:),dx, N,IMAX)
    end do


    !define boundary conditions for the field needed
    unew(iphi,1) = 0.

    !integrate radially filling from point I0 to IMAX
    !since the rhs evaluation might have to be evaluated to
    !high enough order, we have an interpolation, but when
    !trying first or last point there aren't enough points
    !around to do so. Thus we introduce the offset ll which
    !is needed only at those points. 


    !also the stuff below is highly inefficient,
    !i am not discriminating on the integration, as I am integrating
    !all variables. This is horrible, as only one field is being integrated
    !but this can be easily fixed, do not want to worry about it
    I0 = 2

    do i=I0,IMAX
       ll  = 0
       if(i.eq.2.or.i.eq.IMAX) then
          ll = 1
       end if

       del = dx

       uguess = unew(:,i-1)
       dum = du(:,i-1)
       um = unew(:,i-1)
       xm=x(i-1)
       call rhs_space(uguess,um,dum,d2u(:,i-1),k1,del,N,x(i-1),time)
       u1 = unew(:,i-1)+0.5*del*k1

       dum(:) = ( 9./16.*(du(:,i-1)+du(:,i)) - 1./16.* (du(:,i-2+ll)+du(:,i+1-ll)) )
       um(:) = ( 9./16.*(unew(:,i-1)+unew(:,i)) - 1./16.* (unew(:,i-2+ll)+unew(:,i+1-ll)) )
       xm = 0.5*(x(i-1)+x(i))
       call rhs_space(u1,um,dum,d2u(:,i-1),k2,del,N,xm,time)
       u2 = unew(:,i-1)+0.5*del*k2

       dum(:) = ( 9./16.*(du(:,i-1)+du(:,i)) - 1./16.* (du(:,i-2+ll)+du(:,i+1-ll)) )
       um(:) = ( 9./16.*(unew(:,i-1)+unew(:,i)) - 1./16.* (unew(:,i-2+ll)+unew(:,i+1-ll)) )
       xm = 0.5*(x(i-1)+x(i))
       call rhs_space(u2,um,dum,d2u(:,i-1),k3,del,N,xm,time)
       u3 = unew(:,i-1)+del*k3

       uguess = u3
       dum = du(:,i)
       um = unew(:,i)
       xm=x(i)
       call rhs_space(u3,um,dum,d2u(:,i-1),k4,del,N,x(i),time)

       !here is were i am careful what i write over
       unew(iphi,i) = unew(iphi,i-1)+del/6.*(k1(iphi)+2.*k2(iphi)+2.*k3(iphi)+k4(iphi))


    end do

    deallocate(d2u)

  end subroutine rk_fullradial

  
  ! Fix the g boundary condition after doing a time step
  subroutine boundary_fix(u,x,dx,time,dt,ak,tag)
    use m_pars, only : ig,imax

    implicit none
    real*8, dimension(:,:) ::  U
    integer tag
    real*8 time, dt, dx, ak
    real*8, dimension(:) :: x

    !this is here just in case it is needed
    u(ig,IMAX) = u(ig,1)

  end subroutine boundary_fix


end module M_EVOLVE
