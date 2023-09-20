module global
   use, intrinsic :: iso_c_binding
   use constants
   use types

   implicit none

   !  Coil positions
   integer(8) :: num_coils
   real(8), allocatable :: coil_positions(:,:)
   type(vector_grid), allocatable, dimension(:) :: r_coil
   real(8), allocatable, dimension(:,:,:,:) :: abs_r_coil, inv_abs_r_coil_squared

   !  Geometry and magnetic field

   !  Mesh dimensions
   integer(8) :: len_s, len_th, len_ph
   !  Boozer coordinates mesh
   real(8), allocatable, dimension(:) :: s, th, ph
   real(8) :: delta_s, delta_th, delta_ph
   !  Profiles
   real(8), allocatable :: iota(:)
   !  Magnetic field
   type(vector_grid) :: b_super, xyz_grid
   real(8), allocatable, dimension(:,:,:) :: mod_b, sqrt_g, inv_mod_b2, vol_element
   !  Basis and metric tensor
   type(vector_grid) :: e_sub_s, e_sub_th, e_sub_ph
   type(metric_tensor) :: g_sub_ij

   !  Potentials
   !  Time array
   integer(8) :: len_t
   real(8), allocatable, dimension(:) :: t
   !  Potential fourier expansion
   integer(8) :: num_modes
   integer(8), allocatable, dimension(:) :: ms_pot, ns_pot
   real(8), allocatable, dimension(:) :: ws_pot
   complex(8), allocatable, dimension(:,:) :: potnm
   !  Complex exponentials for potential expansion
   complex(8), allocatable, dimension(:,:) :: exp_mth_pot, exp_nph_pot, exp_wt_pot
   !  Derivatives of the potential
   real(8), allocatable, dimension(:,:,:) :: dpot_dth, dpot_dph
   !
   type(vector_grid) :: gradpar_pot_super, gradpar_pot_sub, curl_curl_gradpar_pot

   !  Gradients etc
   !  FFTfreqs for derivatives
   real(8), allocatable, dimension(:) :: th_freqs, ph_freqs
   real(8), allocatable, dimension(:,:):: ftheta, fphi
   !  FFT plans
   type(C_PTR) :: fft_plan_2d, ifft_plan_2d

contains

end module global
