module pde_routines
  use constants, only: pi
  implicit none

contains

!!!!!!!!!!!!!! USER INPUT FUNCTIONS HERE !!!!!!!!!!!!!!!!!!!

  ! Function defining the t=0 initial condition
  function init(x) result(y)
    real :: x,y
    y = sin(2*pi*x)
  end function init

  ! Function defining the initial t=0 initial derivatives
  function d_init(x) result(y)
    real :: x,y
    y = 0.*x!5.*pi*cos(1*pi*x)
  end function d_init

  ! Function defining the potential to be considered
  function potent(x) result(y)
    real :: x,y
    y = 0.*x ! sin(1*pi*x)
  end function potent
!!!!!!!!!!!!!! USER INPUT FUNCTIONS END HERE !!!!!!!!!!!!!!!!
  ! Function for calculating the x position for index
  function x_tran(j,nx,dx) result(x)
    real :: dx, x
    integer :: j, nx

    x = (j-nx/2)*dx

  end function x_tran


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
  subroutine initial1(psi1, psi0, nx, dx, dt)
    ! Declaring the types,
    real, intent(inout) :: psi1(0:nx), psi0(0:nx)
    integer, intent(in) :: nx
    real :: dx, dt, x

    ! Declare halfway variables
    real :: ini1, ini2, ini3

    integer :: j

    real :: a

    a = (dt/dx)**2

    do j=1,nx-1

      ! Calculate the x coordinate

      ! Need a bit more of a complicated process
      ini1 = init(j*dx)+dt*d_init(j*dx)
      ini2 = (1./2)*((dt/dx)**2)*(-init((j+1)*dx) &
              +2*init(j*dx)-init((j-1)*dx))

      ! Note that this term contains the potential
      ini3 = (1./2)*((dt/dx)**2)*(potent(j*dx)*init(j*dx))

      psi1(j) = ini1 + ini2 - ini3
    end do

    ! TODO fix this trash
    psi1(0) = 0.
    psi1(nx) = 0.

  end subroutine initial1

  ! Subroutine that determines the value of a point in the
  ! next timestep
  subroutine scheme(psi0, psi1, psi2, nx, dx, dt, j)
    ! Declaring types
    real, intent(in) :: psi0(0:nx)
    real, intent(in) :: psi1(0:nx)
    real, intent(inout) :: psi2(0:nx)

    integer :: nx, i, j
    real :: dx, dt, a

    ! determine the value of a (alpha)
    a = (dt/dx)**2

    ! Calculate the next value
    psi2(j) = a*(psi1(j+1)+psi1(j-1)) &
    + (2*(1-a)  - potent(j*dx)*dt**2)*psi1(j)- psi0(j)

    ! Set the boundaries TODO move this trash
    psi2(0) = 0.
    psi2(nx) = 0.

  end subroutine scheme

end module pde_routines
