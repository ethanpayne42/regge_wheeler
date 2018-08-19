module inout

  implicit none

contains

  ! Function for reading in the user setup
  subroutine read_input(nx, dx, nt, dt)
    integer, intent(inout) :: nx, nt
    real, intent(inout) :: dx, dt

    ! Ask for the number of points
    ! for the radial direction
    print*,'Input number of points in x direction'
    print*,'and spacing in space'
    read (*,*) nx, dx
    print*,'We have a total spatial length of',nx*dx
    print*,'Input number of points in t direction'
    print*,'and spacing in time'
    read (*,*) nt,dt

  end subroutine read_input

end module inout
