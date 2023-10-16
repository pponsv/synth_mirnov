module constants

   use iso_fortran_env

   implicit none

   !  Kinds
   integer, parameter :: r8 = real64
   integer, parameter :: i8 = int64

   !  Physical constants
   real(r8), parameter :: PI   = 4.0_r8 * atan(1.0_r8)
   real(r8), parameter :: C    = 299792458.0_r8
   real(r8), parameter :: EPS0 = 8.8541878128E-12_r8
   real(r8), parameter :: MU0  = 1.25663706212E-6_r8

   !  Misc
   complex(r8), parameter :: imag = cmplx(0, 1.0_r8, kind=r8)

end module constants
