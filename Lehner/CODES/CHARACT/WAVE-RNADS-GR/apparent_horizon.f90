!This module calculates the apparent horizon on every time slice
!and outputs the rho_ah and lets rho_out be 10 points inside the AH

MODULE M_APPARENT_HORIZON

  implicit none
  
contains
  
  
  subroutine apparent_horizon(s_hat,sigma_hat,x,N_total,rho_ah,rho_out_new,area)

     use M_pars, only: rads, N, mask, maskhor, IMAX, IMAXOLD

     implicit none

     integer :: i
     real*8, intent(out) :: rho_ah, rho_out_new, area
     integer, intent(out) :: N_total
     real*8 :: y1, y2, L, sigma_inter
     real*8, dimension(:), intent(inout) :: s_hat, sigma_hat
     real*8, dimension(:), intent(inout) :: x


     L = rads
     i = 1
     do while (s_hat(i)+(1./(2.0*L**2.))*(sigma_hat(i)&
          &+1/x(i))**2. +0.5 .GT. 0 .AND. i .LE. N)
        i = i+1
     end do

!	IMAX= i


     y1 = s_hat(i)+(1./(2.0*L**2.))*(sigma_hat(i)+1/x(i))**2. + 0.5
     y2 = s_hat(i-1)+(1./(2.0*L**2.))*(sigma_hat(i-1)+1./x(i-1))**2. + 0.5

     rho_ah = x(i) - y1 * ((x(i)-x(i-1)) / (y1-y2) )

     sigma_inter = (rho_ah - x(i))*(((sigma_hat(i))&
          &-(sigma_hat(i-1)))/(x(i)-x(i-1))) &
          &+ (sigma_hat(i))
     area = 4. * 3.14159265358979 * (sigma_inter + 1./rho_ah)**2.

     !Add ten points inside the horizon
!force it not to go backwards just in case
!     if(i + 10.lt.IMAX) 
      if(i+10.gt.IMAXOLD) N_total = min(i+10,IMAXOLD + 1)
      if(i+10.lt.IMAXOLD) N_total = max(i+10,IMAXOLD - 1)
      if(i+10.eq.IMAXOLD) N_total = IMAXOLD

!make sure we do not go too far. 
     N_total = min(N_total,N)

     rho_out_new = x(N_total)

     IMAX = N_total
     mask(N_total:N) = 0.0d0
     maskhor(i:N) = 0.0d0

!commenting this out, i used it to check indeed
!the masked region never goes inwards
!    write(211,*) IMAX

     return

   end subroutine apparent_horizon



   subroutine area_evolution(area, dt, darea)

     integer :: j
     real*8 :: dt
     real*8, dimension(:), intent(inout) :: area, darea


     !first and last point using second order FFD and BFD, CFD middle
     darea(1) = (-1.5*area(1) + 2*area(2) - 0.5*area(3)) / dt

     darea(2:size(area)-1) = (/ ( (area(j+1)-area(j-1)) / (2*dt), j=2,(size(area)-1)) /)

     darea(size(area)) = (0.5*area(size(area)-2) - 2*area(size(area)-1) &
          &+ 1.5*area(size(area))) / dt


     return

   end subroutine area_evolution

END MODULE M_APPARENT_HORIZON
