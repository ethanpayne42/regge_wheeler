module inout

  implicit none

contains

  ! Function for reading in the user setup
  subroutine read_input(nx, lx, nt, lt)
    integer, intent(inout) :: nx, nt
    real, intent(inout) :: lx, lt

    ! Ask for the number of points
    ! for the radial direction
    print*,'Input number of points in x direction'
    print*,'and the total length in M'
    read (*,*) nx,lx
    print*,'Input number of points in t direction'
    print*,'and the total length in time in 1./M'
    read (*,*) nt,lt

  end subroutine read_input

end module inout
