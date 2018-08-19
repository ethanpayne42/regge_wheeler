module pde_routines
  use constants, only: pi
  implicit none

contains

!!!!!!!!!!!!!! USER INPUT FUNCTIONS HERE !!!!!!!!!!!!!!!!!!!

  ! Function defining the t=0 initial condition
  function init(x) result(y)
    real :: x,y
    y = exp(-(5*x)**2)
  end function init

  ! Function defining the initial t=0 initial derivatives
  function d_init(x) result(y)
    real :: x,y
    y = 0.*x
  end function d_init

  ! Function defining the potential to be considered
  function potent(x) result(y)
    real :: x,y
    y = 0.*x
  end function potent
!!!!!!!!!!!!!! USER INPUT FUNCTIONS END HERE !!!!!!!!!!!!!!!!
  ! Function for calculating the x position for index
  function x_tran(j,dx,nx) result(x)
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
    real :: x0

    integer :: j

    do j=0,nx

      ! Calculate the x coordinate
      x0 = x_tran(j, dx, nx)

      ! Function for the initial condition
      psi0(j) = init(x0)

    end do

  end subroutine initial0

  ! Subroutine for the first timestep given the initial
  ! conditions, note that it is the same order as the rest
  subroutine initial1(psi1, nx, dx, dt)
    ! Declaring the types,
    real, intent(inout) :: psi1(0:nx)
    integer, intent(in) :: nx
    real :: dx, dt

    ! Declare halfway variables
    real :: ini1, ini2, ini3

    integer :: j

    real :: a
    real :: x0, x1, xn1

    a = (dt/dx)**2

    do j=1,nx-1

      ! Calculate the x coordinate
      x0 = x_tran(j, dx, nx)
      x1 = x_tran(j+1, dx, nx)
      xn1 = x_tran(j-1, dx, nx)

      ! Need a bit more of a complicated process
      ini1 = init(x0)+dt*d_init(x0)
      ini2 = (1./2)*((dt/dx)**2)*(-init(x1) &
              +2*init(x0)-init(xn1))

      ! Note that this term contains the potential
      ini3 = (1./2)*((dt/dx)**2)*(potent(x0)*init(x0))

      psi1(j) = ini1 + ini2 - ini3
    end do

  end subroutine initial1

  ! Here is the subroutine for the
  ! Absorbing boundary conditions (Mur ABC)
  ! which should stop the wave from
  ! propagating back


  ! Subroutine that determines the value of a point in the
  ! next timestep
  subroutine scheme(psi0, psi1, psi2, nx, dx, dt, j)
    ! Declaring types
    real, intent(in) :: psi0(0:nx)
    real, intent(in) :: psi1(0:nx)
    real, intent(inout) :: psi2(0:nx)

    integer :: nx, j
    real :: dx, dt, a

    real :: x0

    ! Calculate the x coordinate
    x0 = x_tran(j, dx, nx)

    ! determine the value of a (alpha)
    a = (dt/dx)**2

    ! Calculate the next value
    psi2(j) = a*(psi1(j+1)+psi1(j-1)) &
    + (2*(1-a)  - potent(x0)*dt**2)*psi1(j)- psi0(j)

    ! Set the boundaries TODO move this trash
    psi2(0) = 0.
    psi2(nx) = 0.

  end subroutine scheme

end module pde_routines
