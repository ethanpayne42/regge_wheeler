module pde_routines
  use constants, only: pi
  implicit none

contains

!!!!!!!!!!!!!! USER INPUT FUNCTIONS HERE !!!!!!!!!!!!!!!!!!!

  ! Function defining the t=0 initial condition
  function init(x) result(y)
    real :: x,y
    y = sin(5*pi*x)
  end function init

  ! Function defining the initial t=0 initial derivatives
  function d_init(x) result(y)
    real :: x,y
    y = (5*pi)*cos(5*pi*x)
  end function d_init

  ! Function defining the potential to be considered
  function potent(x) result(y)
    real :: x,y
    y = 0*x
  end function potent
!!!!!!!!!!!!!! USER INPUT FUNCTIONS END HERE !!!!!!!!!!!!!!!!

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

  ! Subroutine for the first timestep given the initial
  ! conditions, note that it is the same order as the rest
  subroutine initial1(psi1, nx, dx, dt)
    ! Declaring the types,
    real, intent(inout) :: psi1(0:nx)
    integer, intent(in) :: nx
    real, intent(in) :: dx, dt

    ! Declare halfway variables
    real :: ini1, ini2, ini3

    integer :: j

    do j=0,nx
      ! Need a bit more of a complicated process
      ini1 = init(j*dx)+dt*d_init(j*dx)
      ini2 = (1./2)*((dt/dx)**2)*(init((j+1)*dx) &
              -2*init(j*dx)+init((j-1)*dx))

      ! Note that this term contains the potential
      ini3 = (1./2)*((dt/dx)**2)*(potent(j*dx)*init(j*dx))

      psi1(j) = ini1 + ini2 - ini3

    end do

  end subroutine initial1

end module pde_routines
