module constants

   implicit none

   !  Kinds
   integer, parameter :: r8 = selected_real_kind(15, 200)
   integer, parameter :: i8 = selected_int_kind(16)

   !  Physical constants
   real(r8), parameter :: PI   = 4.0_8 * atan(1.0_8)
   real(r8), parameter :: C    = 299792458.0_8
   real(r8), parameter :: EPS0 = 8.8541878128e-12_8
   real(r8), parameter :: MU0  = 1.25663706212e-6_8

   !  Misc
   complex(r8), parameter :: imag = complex(0, 1._r8)

end module constants
