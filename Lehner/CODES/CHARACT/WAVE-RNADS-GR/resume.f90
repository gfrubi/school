module M_resume
  implicit none
  
contains

  subroutine resume(uold, time, n)
    use m_pars, only : isigma, itau, is, ialpha, ibeta, iw, iz, ipir, ipii, iphir, iphii, iprefix, resumelevel,mask, IMAX
    implicit none
    real*8, dimension(:,:), intent(inout) :: uold
    real*8, intent(out) :: time
    real*8, dimension(:), allocatable :: coordsgft
    integer, dimension(1) :: shapegft
    character(len=32) :: cnamesgft
    integer :: ret, gft_read_full, gft_read_brief, level, i, n

    allocate(coordsgft(n))


    !initialize all to zero, then fill
    uold(:,:) = 0.0d0
    
    ! If resumelevel=-m, then resume from the mth level from the end
    ! of the sdf file.  Cycle through to figure out how many levels
    ! there are.
    if (resumelevel < 0) then
       i=1
       do
          ret = gft_read_full(trim(iprefix)//'phir',i,shapegft,cnamesgft,1,time,coordsgft,uold(iphir,:))
          if (ret < 1) exit
          i=i+1
       end do
       level=i+resumelevel
    else
       level=resumelevel
    end if

    ! Fill up other fields
    ret = gft_read_brief(trim(iprefix)//'phii',level,uold(iphii,:))
    ret = gft_read_brief(trim(iprefix)//'pir',level,uold(ipir,:))
    ret = gft_read_brief(trim(iprefix)//'pii',level,uold(ipii,:))
    ret = gft_read_brief(trim(iprefix)//'sigma',level,uold(isigma,:))
    ret = gft_read_brief(trim(iprefix)//'tau',level,uold(itau,:))
    ret = gft_read_brief(trim(iprefix)//'alpha',level,uold(ialpha,:))
    ret = gft_read_brief(trim(iprefix)//'beta',level,uold(ibeta,:))
    ret = gft_read_brief(trim(iprefix)//'w',level,uold(iw,:))
    ret = gft_read_brief(trim(iprefix)//'z',level,uold(iz,:))
    ret = gft_read_brief(trim(iprefix)//'s',level,uold(is,:))
    ret = gft_read_brief(trim(iprefix)//'mask',level,mask)

    print *, "Loaded data at time =", time

    ! Set IMAX from looking at the mask array
    IMAX = n
    do i=1, n
       if (mask(i).lt.0.01) then
          IMAX = i
          exit
       end if
    end do

    deallocate(coordsgft)

  end subroutine resume

end module M_resume
