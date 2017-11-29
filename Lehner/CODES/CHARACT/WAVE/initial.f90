	MODULE M_INITIAL
        use M_pars
        use M_derivs
	implicit none
	
	CONTAINS
	
	subroutine initial(Uold,dx,x,Nv,N,dt,lambda)
	use m_pars, only : iphi, ipi, ib,ibeta,ia, ialp, disip, amp
        implicit none
	real*8, dimension(:,:) :: Uold
	real*8, dimension(:) :: x
	real*8 :: dx, dy,dt,lambda
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

	do i=1,N
	  if(x(i).ge.6.and.x(i).lt.8) then
	     uold(ig,i) = amp *(x(i)-6.)**3*(x(i)-8.)**3*(2.*x(i)-6.-8.)
	  end if

 	end do


	end subroutine initial
	
	
	END MODULE M_INITIAL