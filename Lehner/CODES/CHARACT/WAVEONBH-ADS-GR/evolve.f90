!this gives and error that goes ~ dt^2 , dx^4, so if 
!the spacing is sufficiently small, one sees the dx^4 error...

module M_EVOLVE
  implicit none

contains

  ! try getting higher order, rk in time and space in 
  ! sequential order 

  subroutine evolve(unew,uold,dx,dt,time,x,Nv,Nspatial,N,lambda)
    use m_derivs
    use m_boundary
    use m_rhs
    use m_pars, only : isigma, itau, is, ialpha, ibeta, ipi, iphi, disip, &
         &     a2iold, a2i, p2, p2old, amp, IMAX, a2, mask, freqhor,initime,&
         &     area_rad, IMAXOLD, null_rad, ITE

    implicit none

    real*8, dimension(:,:), intent(inout):: unew
    real*8, dimension(:,:), intent(in):: uold
    real*8, dimension(:):: x	
    real*8 dx, dt, time,lambda,rback
    integer :: Nv,N,Nspatial

    !here are some basic local variables, modify as needed

    real*8 :: xm, dxPhim, Phim, A2xx, A1x, A1t, A0, Fac
    real*8, dimension(:,:), allocatable :: tempdis, u1,u2,u3,u4,du
    real*8, dimension(:,:), allocatable :: urk1,urk2,urk3,urk4
    real*8 :: ak1, ak2, ak3, ak4, akt

    real*8 EPS, factor2,factor3,p20, coup
    integer i,j,gft_out_full,gft_out_brief, ret
    logical :: ltrace
    real*8 :: coeff(3), p0, dtp0, d2tp0
    real*8 :: tangh, Axx, aa, bb, bet, d2t_p0, dt_p0, s1, AL, AR, XS,BS,AS,ATxx


    coup = 1.0
    !for debug
    ltrace = .false.
    !allocate memory
    allocate(tempdis(Nv,N),u1(Nv,N),u2(Nv,N),u3(Nv,N),u4(Nv,N),du(Nv,N), &
         &           urk1(Nv,N),urk2(Nv,N),urk3(Nv,N),urk4(Nv,N))

    !integrate here so that the fields are field with proper values initially
    if(time.lt.initime+0.25*dt) then
       !call boundary_fix(uold,x,dx,time,dt,ak1,0)
       call rk_fullradial(uold,du,x,dx,dt,Nv,Nspatial,N,0,time)
       ! Output initial data
       ret = gft_out_full('phi',time, n, 'x', 1, x, Uold(iphi,:))
       ret = gft_out_full('pi',time, n, 'x', 1, x, Uold(ipi,:))
       ret = gft_out_full('s',time, n, 'x', 1, x, Uold(is,:))
       ret = gft_out_full('alpha',time, n, 'x', 1, x, Uold(ialpha,:))
       ret = gft_out_full('sigma',time, n, 'x', 1, x, Uold(isigma,:))
       ret = gft_out_full('beta',time, n, 'x', 1, x, Uold(ibeta,:))
       ret = gft_out_full('tau',time, n, 'x', 1, x, Uold(itau,:))
    end if

    !LUIS
    u1 = uold
    u2 = uold
    u3 = uold
    u4 = uold

    !march up 1st step rk, then full radial out
    !notice the rk evoln step is always the same, only changing
    !by the size of the step, which we pass to that routine
    call rk_timeup(uold,uold,u1,du,tempdis,urk1,x,dx,0.5*dt,Nv,Nspatial,N,time)
    !call boundary_fix(u1,x,dx,time,dt*0.5,ak1,0)
    call rk_fullradial(u1,du,x,dx,dt,Nv,Nspatial,N,0,time+0.5*dt)

    !print *, u1(ipi,:)
    !stop
    
    !now do the 2nd RK step, then full radial out
    call rk_timeup(uold,u1,u2,du,tempdis,urk2,x,dx,0.5*dt,Nv,Nspatial,N,time+0.5*dt)
    !call boundary_fix(u2,x,dx,time+0.5*dt,dt*0.5,ak2,0)
    call rk_fullradial(u2,du,x,dx,dt,Nv,Nspatial,N,0,time+0.5*dt)


    !now do the 3rd RK step, then full radial out
    call rk_timeup(uold,u2,u3,du,tempdis,urk3,x,dx,dt,Nv,Nspatial,N,time+dt)
    !call boundary_fix(u3,x,dx,time+dt,dt,ak3,0)
    call rk_fullradial(u3,du,x,dx,dt,Nv,Nspatial,N,0,time+dt)


    !now do the 4th RK step, then full radial ou
    call rk_timeup(uold,u3,u4,du,tempdis,urk4,x,dx,dt,Nv,Nspatial,N,time+dt)
    !call boundary_fix(u4,x,dx,time+dt,dt,ak4,0)

    !now combine the RKs
    unew(iphi,1:IMAX)=uold(iphi,1:IMAX) &
         & + dt/6.*(urk1(iphi,1:IMAX) + 2.*urk2(iphi,1:IMAX) &
         &          + 2.*urk3(iphi,1:IMAX) + urk4(iphi,1:IMAX))
    !unew(ig,IMAX)=unew(ig,1)


    call rk_fullradial(unew,du,x,dx,dt,Nv,Nspatial,N,0,time+dt)


    deallocate(tempdis,u1,u2,u3,u4,du,urk1,urk2,urk3,urk4)


    !FIND HORIZON
    IF(mod(ITE,freqhor).eq.0) THEN
       !here do something special if finding a horizon is required
       !but do so every so often to not waste too much time

       !do something so that it does not move more than 1 pt at a time

    END IF
    !END FIND HORIZON

    !null horizon

!!$    do i=1, N
!!$       xm = x(i)
!!$       !	Axx = xm**2*p0**2/12+xm**4*log(xm)*p0**4/72-xm**2*unew(ia,i)/2+xm**4*log(xm&
!!$       !     &)*p0*d2tp0/12-xm**4*log(xm)*dtp0**2/12-1.D0/2.D0
!!$       !use horizon definition,i'll still call it Axx being lazy 
!!$       aa = unew(ia,i)
!!$       bb = unew(ib,i)
!!$       bet = unew(ibeta,i)
!!$
!!$       Axx = 1-xm**2*p0**2/6+xm**4*log(xm)*dt_p0**2/6-xm**4*log(xm)*p0*d2t&
!!$            &_p0/6-xm**4*log(xm)*p0**4/36+xm**2*aa
!!$
!!$       AR = Axx
!!$       if(Axx.lt.0) then
!!$          !          print*,xm, 'hor'
!!$          !use a linear interp to get the right value
!!$          if(abs(AL-AR).gt.1e-8) then
!!$             XS = (x(i)*AL-x(i-1)*AR)/(AL-AR)
!!$          else
!!$             XS = x(i)
!!$          end if
!!$          BS= (unew(ib,i-1)*(XS-x(i))-unew(ib,i)*(XS-x(i-1)))/dx
!!$          null_rad = -1./12.*XS**2*p0**2-1./9.*XS**3*p0*dtp0+XS**2*BS
!!$          null_rad = exp(area_rad)/XS
!!$          exit
!!$       end if
!!$       AL=AR
!!$    end do

    return

  end subroutine evolve

  ! RK intermidiate step for the time advance
  subroutine rk_timeup(uold,uprev,unew,du,tempdis,urk,x,dx,dt,Nv,Nspatial,N,time)
    use m_derivs
    use m_boundary
    use m_rhs
    use m_pars, only : iphi, disip, IMAX, mask
    implicit none
    real*8, dimension(:,:) :: uold, uprev,unew, du,urk,tempdis
    real*8, dimension(:) :: x
    real*8 :: dx, dt, time
    integer :: Nv,Nspatial,N

    integer :: i, j
    real*8 :: coup
    real*8 :: rback

    coup = 1.

    ! We need dissipation for time-evolved fields
    do j=Nspatial+1,Nv
       call dissip(uprev(j,:),tempdis(j,:),dx, N,1,IMAX)
    end do

    ! Equations involve spatial derivatives of phi only
    call derivs(uprev(iphi,:),du(iphi,:),dx, N,1,IMAX)

    
    call rhs_time(uprev,du,urk,dx,N,x,time)

    tempdis(iphi,1) = 0.0d0

    urk = urk + sign(1.0,dx)*disip*tempdis

    do i = 1, IMAX         
       unew(iphi,i) = uold(iphi,i) + dt * urk(iphi,i) 
    end do

  end subroutine rk_timeup

  !RK full march radially
  subroutine rk_fullradial(unew,du,x,dx,dt,Nv,Nspatial,N,pp,time)
    use m_derivs
    use m_boundary
    use m_rhs
    use m_pars, only : isigma, itau, is, ialpha, ibeta, ipi, iphi, disip, &
         &                    a2iold, a2i, p2, p2old, amp,lambda, IMAX, a2, madm
    implicit none
    real*8, dimension(:,:) :: unew, du
    real*8, dimension(:) :: x
    real*8 :: dx, dt,time
    integer :: Nv, Nspatial, N, pp, ll 

    real*8 :: xm,del
    real*8, dimension(Nv) :: k1,k2,k3,k4, u1,u2,u3,u4, dum,uguess, um
    integer :: i,j,I0

    integer :: gft_out_full,gft_out_brief, ret, coun, ic, INFO
    real*8, dimension(:,:), allocatable :: d2u
    real*8 :: coeff(3), p0, dtp0, d2tp0, tangh, pixx

    allocate(d2u(Nv,N))


    !below is a call to derivatives that would be needed
    !if the rhs of the radial equation has them 
    !in this toy model they are not needed, 
    !but will leave here in case they are needed
    !in the future      
!!$    do j=1,Nv
!!$       call derivs(unew(j,:),du(j,:),dx, N,1,IMAX)
!!$       call dersecond2(unew(j,:),d2u(j,:),dx, N,IMAX)
!!$    end do

    ! We only need first derivatives of phi field. Set all others to
    ! 0, but we won't use them.
    du = 0.
    d2u = 0.
    call derivs(unew(iphi,:),du(iphi,:),dx, N,1,IMAX)

    ! Outer boundary conditions for the fields
    ! The beta field takes the ADM mass as the boundary condition
    unew(:Nspatial,1) = 0.
    unew(ibeta,1) = -2*madm

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
       call rhs_space(uguess,um,dum,d2u(:,i-1),k1,del,Nv,Nspatial,N,x(i-1),time)
       u1 = unew(:,i-1)+0.5*del*k1

       dum(:) = ( 9./16.*(du(:,i-1)+du(:,i)) - 1./16.* (du(:,i-2+ll)+du(:,i+1-ll)) )
       um(:) = ( 9./16.*(unew(:,i-1)+unew(:,i)) - 1./16.* (unew(:,i-2+ll)+unew(:,i+1-ll)) )
       xm = 0.5*(x(i-1)+x(i))
       call rhs_space(u1,um,dum,d2u(:,i-1),k2,del,Nv,Nspatial,N,xm,time)
       u2 = unew(:,i-1)+0.5*del*k2

       dum(:) = ( 9./16.*(du(:,i-1)+du(:,i)) - 1./16.* (du(:,i-2+ll)+du(:,i+1-ll)) )
       um(:) = ( 9./16.*(unew(:,i-1)+unew(:,i)) - 1./16.* (unew(:,i-2+ll)+unew(:,i+1-ll)) )
       xm = 0.5*(x(i-1)+x(i))
       call rhs_space(u2,um,dum,d2u(:,i-1),k3,del,Nv,Nspatial,N,xm,time)
       u3 = unew(:,i-1)+del*k3

       uguess = u3
       dum = du(:,i)
       um = unew(:,i)
       xm=x(i)
       call rhs_space(u3,um,dum,d2u(:,i-1),k4,del,Nv,Nspatial,N,x(i),time)

       !here is were i am careful what i write over
       unew(:Nspatial,i) = unew(:Nspatial,i-1)+del/6.*(k1(:Nspatial)+2.*k2(:Nspatial)+2.*k3(:Nspatial)+k4(:Nspatial))
       
    end do

    !print *, "Output from rk_fullradial"
    !print *,unew(ipi,:)


    !        if(pp.eq.1) then
    !          ret = gft_out_brief('iphi',0.0, (/n/), 1, unew(iphi,:))  
    !	STOP
    !         end if

    deallocate(d2u)

  end subroutine rk_fullradial



  !this subroutine gets a2 and related

  subroutine boundary_fix(u,x,dx,time,dt,ak,tag)
    use m_pars, only : iphi,imax,a2,a2i,a2iold,p2,p2old, amp, lambda,freq

    implicit none
    real*8, dimension(:,:) ::  U
    integer N, Nv, tag
    real*8 time, dt, dx, ak
    real*8, dimension(:) :: x

    real*8 :: coeff(3), p0, dtp0, d2tp0
    real*8 :: tangh

    !this is here just in case it is needed
    !u(ig,IMAX) = u(ig,1)

  end subroutine boundary_fix


  !this subroutine fits using 20 pts
  !also in case it is needed
  subroutine fit(yin,xin,coeff,dx)
    implicit none
    real*8,dimension(:) :: yin, xin
    real*8 :: y(20), x(20), coeff(3)
    real*8 :: dx
    integer :: i

    y = yin(4:23)
    x = xin(4:23)

    coeff(1) = 571.D0/1540.D0*y(1)+459.D0/1540.D0*y(2)+51.D0/220.D0*y(3&
         &)+53.D0/308.D0*y(4)+183.D0/1540.D0*y(5)+111.D0/1540.D0*y(6)+7.D0/2&
         &20.D0*y(7)-3.D0/1540.D0*y(8)-9.D0/308.D0*y(9)-y(10)/20-9.D0/140.D0&
         &*y(11)-111.D0/1540.D0*y(12)-113.D0/1540.D0*y(13)-3.D0/44.D0*y(14)-&
         &87.D0/1540.D0*y(15)-59.D0/1540.D0*y(16)-3.D0/220.D0*y(17)+27.D0/15&
         &40.D0*y(18)+17.D0/308.D0*y(19)+153.D0/1540.D0*y(20)
    coeff(2) = -(6669*y(1)+4827*y(2)+3175*y(3)+1713*y(4)+441*y(5)-641*y&
         &(6)-1533*y(7)-2235*y(8)-2747*y(9)-3069*y(10)-3201*y(11)-3143*y(12)&
         &-2895*y(13)-2457*y(14)-1829*y(15)-1011*y(16)-3*y(17)+1195*y(18)+25&
         &83*y(19)+4161*y(20))/dx/87780
    coeff(3) = (57*y(1)+39*y(2)+23*y(3)+9*y(4)-3*y(5)-13*y(6)-21*y(7)-2&
         &7*y(8)-31*y(9)-33*y(10)-33*y(11)-31*y(12)-27*y(13)-21*y(14)-13*y(1&
         &5)-3*y(16)+9*y(17)+23*y(18)+39*y(19)+57*y(20))/dx**2/17556

  end subroutine fit


end module M_EVOLVE
