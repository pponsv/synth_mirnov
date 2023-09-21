module synthetic_mirnov
   use constants
   implicit none

contains

   subroutine init_booz(s_b, th_b, ph_b, b_mod_b, sqrt_g_b, phi_b_g, iota_b, x, y, z, &
      len_s_b, len_th_b, len_ph_b)
      use global
      use fft_mod, only: fftfreqs, meshgrid, plan_ffts
      integer(8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(8), intent(in) :: s_b(len_s_b), th_b(len_th_b), ph_b(len_ph_b), phi_b_g, iota_b(len_s_b)
      real(8), intent(in), dimension(len_s_b, len_th_b, len_ph_b) :: b_mod_b, sqrt_g_b, x, y, z
      integer :: i

      len_s = len_s_b
      len_th = len_th_b
      len_ph = len_ph_b

      s = s_b
      th = th_b
      ph = ph_b

      delta_s = s(2)-s(1)
      delta_th = th(2) - th(1)
      delta_ph = ph(2) - ph(1)

      th_freqs = 2*PI*fftfreqs(len_th, delta_th)
      ph_freqs = 2*PI*fftfreqs(len_ph, delta_ph)
      allocate(ftheta(len_th, len_ph), fphi(len_th, len_ph))
      call meshgrid(th_freqs, ph_freqs, ftheta, fphi)

      call plan_ffts()

      iota = iota_b
      mod_b = b_mod_b
      inv_mod_b2 = b_mod_b**(-2)
      sqrt_g = sqrt_g_b
      vol_element = sqrt_g * delta_s * delta_th * delta_ph

      allocate(B_super%u1, B_super%u2, B_super%u3, mold=b_mod_b)
      B_super%u1 = 0
      B_super%u3 = 1/(2*pi*sqrt_g)*phi_b_g
      do i=1, len_s
         B_super%u2(i, :, :) = B_super%u3(i, :, :) * iota(i)
      end do

      xyz_grid%u1 = x
      xyz_grid%u2 = y
      xyz_grid%u3 = z

      call gradpar_pot_super%alloc(len_s, len_ph, len_th)

   end subroutine init_booz

   subroutine init_basis(e_sub_s_b, e_sub_th_b, e_sub_ph_b, len_s_b, len_th_b, len_ph_b)
      use global, only : e_sub_s, e_sub_th, e_sub_ph, g_sub_ij
      integer(8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(8), intent(in), dimension(3, len_s_b, len_th_b, len_ph_b) :: e_sub_s_b, e_sub_th_b, e_sub_ph_b

      call e_sub_s%init(e_sub_s_b)
      call e_sub_th%init(e_sub_th_b)
      call e_sub_ph%init(e_sub_ph_b)

      call g_sub_ij%init(e_sub_s, e_sub_th, e_sub_ph)

   end subroutine init_basis

   subroutine init_pot(pot_profiles, ms, ns, fs, time, num_modes, len_time, len_s)
      use potential, only : initialize_potential
      integer(8), intent(in) :: num_modes, len_time, len_s
      integer(8) :: ms(num_modes), ns(num_modes)
      real(8) :: fs(num_modes), time(len_time)
      complex(8) :: pot_profiles(num_modes, len_s)

      call initialize_potential(pot_profiles, ms, ns, fs, time)

   end subroutine init_pot

   subroutine init_coils(xyzs, ncoils)
      use global, only : coil_positions, num_coils, r_coil, xyz_grid, abs_r_coil, &
         inv_abs_r_coil_squared, len_s, len_th, len_ph
      use types, only : dot_vector_grid
      integer(8), intent(in) :: ncoils
      real(8), intent(in) :: xyzs(3, ncoils)
      integer :: idx_coil

      coil_positions = xyzs
      num_coils = ncoils
      allocate(inv_abs_r_coil_squared(num_coils, len_s, len_th, len_ph))
      allocate(r_coil(num_coils))

      !$OMP PARALLEL DO PRIVATE(idx_coil)
      do idx_coil = 1, num_coils
         r_coil(idx_coil)%u1 = xyz_grid%u1 - coil_positions(1,idx_coil)
         r_coil(idx_coil)%u2 = xyz_grid%u2 - coil_positions(2,idx_coil)
         r_coil(idx_coil)%u3 = xyz_grid%u3 - coil_positions(3,idx_coil)
         inv_abs_r_coil_squared(idx_coil, :, :, :) = dot_vector_grid(r_coil(idx_coil), r_coil(idx_coil))
      end do
      !$OMP END PARALLEL DO
      abs_r_coil = sqrt(inv_abs_r_coil_squared)
      inv_abs_r_coil_squared = 1./inv_abs_r_coil_squared
      print *, (sum(abs_r_coil))

   end subroutine init_coils

   subroutine run(num_coils_tmp, len_t_tmp, db_coils)
      use main
      integer(8), intent(in) :: num_coils_tmp, len_t_tmp
      real(8), intent(out) :: db_coils(num_coils_tmp, 3, len_t_tmp)

      call main_loop(db_coils)

   end subroutine run

   function test_fft_main(in) result(out)
      use fft_mod
      real(8) :: in(:,:), nin(size(in,1), size(in, 2))
      complex(8) :: out(size(in,1), size(in, 2)), cout(size(in,1), size(in, 2))

      call plan_ffts()

      ! print *, in
      out = fft_2d(in)
      cout = out
      nin = ifft_2d(cout)

      ! print *, in - nin

   end function test_fft_main

   function get_potnm(ls, nm) result(out)
      use global, only: potnm
      integer :: ls, nm
      real(8) :: out(nm, ls)

      out = potnm
   end function get_potnm

   subroutine test_meshgrid
      use fft_mod, only: meshgrid
      real(8) :: a(2), b(3)
      real(8) :: aa(2,3), bb(2,3)

      a = (/ 1., 2./)
      b = (/ 3., 2., 1./)
      call meshgrid(a, b, aa, bb)

      print *, 'AA'
      print *, aa
      print *, 'BB'
      print *, bb(1,:)
   end subroutine test_meshgrid


end module synthetic_mirnov
