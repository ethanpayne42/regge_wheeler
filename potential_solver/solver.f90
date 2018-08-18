program solver

  ! pde_routines contains all the information
  ! about the scheme and initial conditions
  use pde_routines
  use inout , only: read_input
  implicit none

  ! Variables for the number of points in x
  ! as well as the total length
  integer :: nx, nt
  real :: lx, lt

  ! Lets declare the i-1 and i indices
  real, allocatable :: psi0(:)
  real, allocatable :: psi1(:)

  ! Lets declare the indices too
  integer :: i ! time index
  integer :: j ! position index


  ! Read in the values for the total lengths
  ! and number of points for arrays
  call read_input(nx,lx,nt,lt)

  ! Lets now set the size of psi0 and psi1
  allocate(psi0(0:nx))
  allocate(psi1(0:nt))
  ! Lets force these to be zero to avoid the random
  ! jumps in the allocated memory values
  psi0 = 0.
  psi1 = 0.

  ! Now we need to set up the initial conditions


end program solver
