	MODULE M_constraints
	implicit none
	
	contains

	subroutine evalconstraints(c1,c2,un,uo,uom1,time,dt,dx,Nx,x,Nv)
	use m_derivs
	use m_pars, only : ia, isigma, iphi, lambda

	implicit none
	real*8, dimension(:,:), intent(in) :: un,uo,uom1
	real*8, dimension(:), intent(inout) :: c1,c2
	real*8, dimension(:), intent(in) :: x
	real*8 :: dx, dt, time
	integer :: Nx, Nv
	logical :: ltrace

!local vars
	integer :: i, j
	real*8, dimension(:,:), allocatable :: dxu, dtu
	real*8, dimension(:,:), allocatable :: d2xu, dtdxu, d2tu
        real*8 :: p0,dtp0,d2tp0,d3tp0, tangh, xp


        real*8 :: phi,sigma,aa
        real*8 :: d1phi,d1sigma,d1aa
        real*8 :: d2phi,d2sigma,d2aa
        real*8 :: dtphi,dtsigma,dtaa
        real*8 :: d2tphi,d2tsigma,d2taa
        real*8 :: dtxphi,dtxsigma,dtxaa
        real*8 :: s1,s2,s3,t0

	allocate(dxu(Nv,Nx),d2xu(Nv,Nx),dtu(Nv,Nx),dtdxu(Nv,Nx),d2tu(Nv,Nx))


!define proper sources
        tangh = tanh((time-0.5*dt)/lambda)
        p0 = 0.5+0.5*tangh
        dtp0 = 0.5*(1.-(tangh)**2)/lambda
        d2tp0 = - tangh*(1-tangh**2)/lambda**2
        d3tp0 = (-1.+4.*tangh**2-3.*tangh**4)/lambda**3

	deallocate(dxu,d2xu,dtu,dtdxu,d2tu)

	end subroutine  evalconstraints


	end module M_constraints
