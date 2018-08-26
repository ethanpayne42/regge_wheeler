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

  subroutine read_in_potential(pot, xs, ind)
    ! ind is the parameter for which potential
    integer :: j, ind
    real :: pot(0:)
    real :: xs(0:)
    real :: tmp(5)
    open(unit = 100, file = 'potentials_data.dat', &
         action = 'read')

    do j=0,3999
      read(100,*) tmp
      xs(j) = tmp(1)
      pot(j) = tmp(ind+1)
    end do

  end subroutine read_in_potential

end module inout
