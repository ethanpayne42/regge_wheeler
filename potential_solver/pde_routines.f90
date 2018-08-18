module pde_routines
  use constants, only: pi
  implicit none

contains

  ! Function defining the t=0 initial condition
  function init(x) result(y)
    real :: x,y

    y = sin(5*pi*x)
  end function init

  ! Subroutine for the initial conditions of the wavepacket
  subroutine initial0(psi0,nx,dx)
    ! Declare the intentions
    real, intent(inout) :: psi0(0:nx)
    integer, intent(in) :: nx
    real, intent(in) :: dx

    integer :: j

    do j=0,nx
      ! Function for the initial condition
      psi0(j) = init(j*dx)

    end do

  end subroutine initial0

end module pde_routines
