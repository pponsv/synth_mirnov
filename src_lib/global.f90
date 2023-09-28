module global
   use, intrinsic :: iso_c_binding
   use constants
   use types

   implicit none

   !  Coil positions
   integer(i8) :: num_coils
   real(r8), allocatable :: coil_xyz(:,:)
   type(vector_grid), allocatable, dimension(:) :: e_sub_s_x_r, &
      e_sub_th_x_r, e_sub_ph_x_r

   !  GEOMETRY AND MAGNETIC FIELD

   !  Mesh dimensions
   integer(i8) :: len_s, len_th, len_ph
   !  Boozer coordinates mesh
   real(r8), allocatable, dimension(:) :: s, th, ph
   real(r8) :: delta_s, delta_th, delta_ph
   !  Profiles
   real(r8), allocatable :: iota(:)
   !  Magnetic field
   type(vector_grid) :: b_super, xyz_grid
   !  Scalar quantities
   real(r8), allocatable, dimension(:,:,:) :: mod_b, sqrt_g, inv_sqrt_g, inv_mod_b2
   !  Basis and metric tensor
   type(vector_grid) :: e_sub_s, e_sub_th, e_sub_ph
   type(metric_tensor) :: g_sub_ij

   !  POTENTIALS

   !  Time array
   integer(i8) :: len_t
   real(r8), allocatable, dimension(:) :: t
   !  Potential fourier expansion
   integer(i8) :: num_modes
   integer(i8), allocatable, dimension(:) :: ms_pot, ns_pot
   real(r8), allocatable, dimension(:) :: ws_pot
   complex(r8), allocatable, dimension(:,:) :: prof_nm
   !  Derivatives of the potential
   complex(r8), allocatable, dimension(:,:,:,:) :: dpot_dth, dpot_dph
   !
   type(vector_grid), allocatable, dimension(:) :: gradpar_pot_super, j_super

   !  Gradients etc
   !  FFTfreqs for derivatives
   real(r8), allocatable, dimension(:) :: th_freqs, ph_freqs
   real(r8), allocatable, dimension(:,:):: fth, fph
   !  FFT plans
   type(C_PTR) :: fft_plan_2d, ifft_plan_2d

contains

end module global
