! File containing the constants module
! Simply (for now) includes the value of pi
module constants
  implicit none
  real, parameter :: pi = 4.*atan(1.)
  complex, parameter :: ic = (0.0, 1.0)
end module constants
