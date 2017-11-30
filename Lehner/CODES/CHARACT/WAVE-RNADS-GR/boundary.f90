	MODULE M_BOUNDARY
	
	implicit none
	
	contains
	
	
	subroutine boundary(unew, time,dx,dt,lambda)
	use m_pars, only:  N, gamma, &
	&                  iphir, iphii
	implicit none	
	real*8, dimension(:,:), intent(inout) :: unew
	real*8 dx,dt, time, factor,lambda
	integer i,j

 !set boundary values
 unew(iphir,1) = 0.0
 unew(iphii,1) = 0.0

	end subroutine boundary
	
	
	end module M_BOUNDARY
	
	
	
