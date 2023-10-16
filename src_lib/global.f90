module global

   use constants

   implicit none

   !  Coil positions
   integer(i8) :: num_coils
   real(r8), allocatable :: coil_xyz(:,:)
   real(r8), allocatable, dimension(:,:,:,:) :: us, uth, uph, r_coil

   !  GEOMETRY AND MAGNETIC FIELD

   !  Mesh dimensions
   integer(i8) :: len_s, len_th, len_ph
   !  Boozer coordinates mesh
   real(r8), allocatable, dimension(:) :: s, th, ph
   real(r8) :: delta_s, delta_th, delta_ph
   !  Profiles
   real(r8), allocatable :: iota(:)
   !  Magnetic field
   real(r8), allocatable, dimension(:,:,:,:) :: xyz_grid
   complex(r8), allocatable, dimension(:,:,:,:) :: b_super
   !  Scalar quantities
   real(r8), allocatable, dimension(:,:,:) :: mod_b, sqrt_g, inv_sqrt_g, inv_mod_b2, r_3_sqrtg
   !  Basis and metric tensor
   real(r8), allocatable, dimension(:,:,:,:) :: e_sub_s, e_sub_th, e_sub_ph
   real(r8), allocatable, dimension(:,:,:,:,:) :: g_sub_ij

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
   !  Parallel gradient, currents
   complex(r8), allocatable, dimension(:,:,:,:,:) :: gradpar_pot_super, j_super

   !  GRADIENTS

   !  FFTfreqs for derivatives
   real(r8), allocatable, dimension(:) :: th_freqs, ph_freqs
   real(r8), allocatable, dimension(:,:) :: fth, fph

end module global
