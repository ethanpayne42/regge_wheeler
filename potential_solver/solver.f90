program solver

  ! pde_routines contains all the information
  ! about the scheme and initial conditions
  use pde_routines, only: initial0, initial1, scheme
  use inout , only: read_input
  implicit none

  ! Variables for the number of points in x
  ! as well as the total length
  integer :: nx, nt
  real :: lx, lt
  real :: dx, dt

  ! Lets declare the i-1 and i indices
  real, allocatable :: psi0(:)
  real, allocatable :: psi1(:)
  real, allocatable :: psi2(:)

  ! Lets declare the indices too
  integer :: i ! time index
  integer :: j ! position index

  ! Read in the values for the total lengths
  ! and number of points for arrays
  call read_input(nx,dx,nt,dt)
  lx = nx*dx
  lt = nt*dt

  ! Lets now set the size of psi0 and psi1
  allocate(psi0(0:nx))
  allocate(psi1(0:nx))
  allocate(psi2(0:nx))
  ! Lets force these to be zero to avoid the random
  ! jumps in the allocated memory values
  psi0 = 0.
  psi1 = 0.
  psi2 = 0.

  ! Now we need to set up the initial conditions
  call initial0(psi0, nx, dx)
  call initial1(psi1, psi0, nx, dx, dt)

  !print*,psi0
  !print*,psi1

  ! Write the first two arrays
  open(1,file='output.dat')
  write(1,*),psi0
  write(1,*),psi1

  ! Now we are at a point to implement the full scheme, almost
  do i=2,nt
    do j=1, nx-1
      ! Calculate the next timestep
      call scheme(psi0, psi1, psi2, nx, dx, dt, j)
    end do
    ! Write to file
    write(1,*),psi2

    ! Make the current step the new step and the old step
    ! the current step
    psi0 = psi1
    psi1 = psi2

  end do

  close(1)

end program solver
