program solver

  ! pde_routines contains all the information
  ! about the scheme and initial conditions
  use pde_routines, only: initial0, initial1, scheme, &
                          x_tran, boundaries
  use inout , only: read_input
  implicit none

  ! Variables for the number of points in x
  ! as well as the total length
  integer :: nx, nt
  real :: lx, lt
  real :: dx, dt

  ! Lets declare the i-1 and i indices
  complex, allocatable :: psi0(:)
  complex, allocatable :: psi1(:)
  complex, allocatable :: psi2(:)

  ! Lets declare the indices too
  integer :: i ! time index
  integer :: j ! position index

  real, allocatable :: xs(:), ts(:)

  ! Read in the values for the total lengths
  ! and number of points for arrays
  !call read_input(nx,dx,nt,dt)

  nx = 100; dx = 2e-2 ! 7000
  nt = 200; dt = 4e-3 ! 50000

  lx = nx*dx
  lt = nt*dt

  print*,'Length in space of', lx
  print*,'Length in time of', lt

  ! Now, we can set up the x and t arrays
  allocate(xs(0:nx))
  allocate(ts(0:nt))

  do j=0,nx
    xs(j) = x_tran(j,dx,nx)
  end do

  ts(0:nt) = (/(i*dt, i=0,nt)/)

  open(100, file='coord.dat')
  write(100,*),xs
  write(100,*),ts
  close(100)

  ! Lets now set the size of psi0 and psi1
  allocate(psi0(0:nx))
  allocate(psi1(0:nx))
  allocate(psi2(0:nx))
  ! Lets force these to be zero to avoid the random
  ! jumps in the allocated memory values
  psi0 = (0.,0.)
  psi1 = (0.,0.)
  psi2 = (0.,0.)

  ! Now we need to set up the initial conditions
  call initial0(psi0, nx, dx)
  call initial1(psi1, nx, dx, dt)

  !print*,psi0
  !print*,psi1

  ! Write the first two arrays
  open(1,file='output_real.dat')
  open(2,file='output_imag.dat')

  ! First lets write the x values across the top
  !write(1,*),0,xs

  write(1,*),real(psi0); write(2,*),aimag(psi0)
  write(1,*),real(psi1); write(2,*),aimag(psi1)

  ! Now we are at a point to implement the full scheme, almost
  do i=2,nt
    do j=1, nx-1
      ! Calculate the next timestep
      call scheme(psi0, psi1, psi2, nx, dx, dt, j)

      call boundaries(psi1,psi2,dx,dt,nx)

    end do

    if (mod(i,100) == 0) then
      write(*,"(A,I6,A,I6)",advance="no") '\b\b\b\b\b\b\b&
      \b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b &
      At step ',i,' of ',nt

    end if

    ! Write to file
    write(1,*),real(psi2); write(2,*),aimag(psi2)


    !print*,psi2

    ! Make the current step the new step and the old step
    ! the current step
    psi0 = psi1
    psi1 = psi2

  end do

  close(1); close(2)

end program solver
