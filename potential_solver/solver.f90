program solver

  ! pde_routines contains all the information
  ! about the scheme and initial conditions
  use pde_routines, only: initial0, initial1, scheme, &
                          x_tran, boundaries, potent
  use inout , only: read_input, read_in_potential
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

  real, allocatable :: ts(:)
  real :: pot(0:4000)
  real :: xs(0:4000)

  integer :: potential_index = 1

  call read_in_potential(pot,xs,potential_index)
  print*, 'Potential at near last position is', pot(3999)

  ! Read in the values for the total lengths
  ! and number of points for arrays
  !call read_input(nx,dx,nt,dt)

  nx = 3999; dx = 1e-1 ! 7000
  nt = 5000; dt = 2.5e-2 ! 45000

  lx = nx*dx
  lt = nt*dt

  print*,'Length in space of', lx
  print*,'Length in time of', lt

  ! Now, we can set up the t array
  allocate(ts(0:nt))

  ! Create the x vector array as well as
  ! save the potential to a file
  open(3,file='potential.dat')

  do j=0,nx
    xs(j) = x_tran(j,dx,nx)
    write(3,*) xs(j), real(potent(xs(j)))
    ! TODO replace the real with complex later
  end do

  print*,'written potential to file'
  close(3)

  ! Construct time array
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

    ! Print the step to the screen
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
