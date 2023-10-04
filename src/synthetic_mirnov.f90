module kinds_f2py

   implicit none

   !  Kinds
   integer, parameter :: r8 = selected_real_kind(15, 200)
   integer, parameter :: i8 = selected_int_kind(16)

end module kinds_f2py

module synthetic_mirnov
   use kinds_f2py
   implicit none

contains

   subroutine init_booz(s_b, th_b, ph_b, b_mod_b, sqrt_g_b, phi_b_g, iota_b, x, y, z, &
      len_s_b, len_th_b, len_ph_b)
      use global
      use fft_mod, only: fftfreqs, meshgrid, plan_ffts
      integer(i8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(r8), intent(in) :: s_b(len_s_b), th_b(len_th_b), ph_b(len_ph_b), phi_b_g, iota_b(len_s_b)
      real(r8), intent(in), dimension(len_s_b, len_th_b, len_ph_b) :: b_mod_b, sqrt_g_b, x, y, z
      integer :: i

      !  Coordinate grids
      len_s = len_s_b
      len_th = len_th_b
      len_ph = len_ph_b

      s = s_b
      th = th_b
      ph = ph_b

      !  We assume a uniform grid. For non-uniform grids, considerable
      !  modification is needed
      delta_s = s(2) - s(1)
      delta_th = th(2) - th(1)
      delta_ph = ph(2) - ph(1)

      !  FFTs for derivatives
      th_freqs = 2*PI*fftfreqs(len_th, delta_th)
      ph_freqs = 2*PI*fftfreqs(len_ph, delta_ph)
      allocate(fth(len_th, len_ph), fph(len_th, len_ph))
      call meshgrid(ph_freqs, th_freqs, fph, fth)
      call plan_ffts

      !  Profiles, magnetic field and jacobian
      iota = iota_b
      mod_b = b_mod_b
      inv_mod_b2 = 1. / b_mod_b**2
      sqrt_g = sqrt_g_b
      inv_sqrt_g = 1. / sqrt_g_b

      allocate(b_super(len_s, len_th, len_ph, 3))
      allocate(xyz_grid(len_s, len_th, len_ph, 3))

      !  phi_b_g is the derivative of the flux wrt the radial coordinate,
      !  equal to the flux at the LCFS if the radial coordinate is s
      b_super(:,:,:,1) = 0
      b_super(:,:,:,3) = -phi_b_g/(2*pi*sqrt_g)
      do i=1, len_s
         b_super(i, :, :, 2) = b_super(i, :, :, 3) * iota(i)
      end do

      !  Cartesian grid
      xyz_grid(:,:,:,1) = x
      xyz_grid(:,:,:,2) = y
      xyz_grid(:,:,:,3) = z

   end subroutine init_booz


   subroutine init_basis(e_sub_s_b, e_sub_th_b, e_sub_ph_b, len_s_b, len_th_b, len_ph_b)
      use global, only : e_sub_s, e_sub_th, e_sub_ph, g_sub_ij, b_super
      use helper, only : metric_tensor, checknans
      integer(i8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(r8), intent(in), dimension(len_s_b, len_th_b, len_ph_b, 3) :: e_sub_s_b, e_sub_th_b, e_sub_ph_b
      integer :: i, j, k

      e_sub_s = e_sub_s_b
      e_sub_th = e_sub_th_b
      e_sub_ph = e_sub_ph_b

      g_sub_ij = metric_tensor(e_sub_s, e_sub_th, e_sub_ph)

   end subroutine init_basis


   subroutine init_pot(profiles, ms, ns, fs, time, num_modes, len_time, len_s)
      use potential, only : init_pot_main => init_pot
      integer(i8), intent(in) :: num_modes, len_time, len_s
      integer(i8), intent(in) :: ms(num_modes), ns(num_modes)
      real(r8), intent(in) :: fs(num_modes), time(len_time)
      complex(r8), intent(in) :: profiles(num_modes, len_s)

      call init_pot_main(profiles=profiles, ms=ms, ns=ns, fs=fs, time=time)

   end subroutine init_pot


   subroutine init_coils(xyzs, ncoils)
      use global, only : coil_xyz, num_coils
      integer(i8), intent(in) :: ncoils
      real(r8), intent(in) :: xyzs(ncoils, 3)

      coil_xyz = transpose(xyzs)
      num_coils = ncoils

   end subroutine init_coils


   subroutine run(num_coils_tmp, len_t_tmp, db_coils)
      use main
      integer(i8), intent(in) :: num_coils_tmp, len_t_tmp
      real(r8), intent(out) :: db_coils(3, num_coils_tmp, len_t_tmp)

      call main_loop(db_coils)

   end subroutine run


   function test_fft_main(in) result(out)
      use fft_mod
      complex(r8) :: in(:,:), nin(size(in,1), size(in, 2))
      complex(r8) :: out(size(in,1), size(in, 2)), cout(size(in,1), size(in, 2))

      call plan_ffts()

      out = fft_2d(in)
      cout = out
      nin = ifft_2d(cout)

   end function test_fft_main


   subroutine test_magnetic_field(i, j, k)
      use global
      integer, intent(in) :: i, j, k

      print *,  b_super(i, j, k, 2) * e_sub_th(i, j, k, 1) + b_super(i, j, k, 3) * e_sub_ph(i, j, k, 1)
      print *,  b_super(i, j, k, 2) * e_sub_th(i, j, k, 2) + b_super(i, j, k, 3) * e_sub_ph(i, j, k, 2)
      print *,  b_super(i, j, k, 2) * e_sub_th(i, j, k, 3) + b_super(i, j, k, 3) * e_sub_ph(i, j, k, 3)

   end subroutine test_magnetic_field


   function get_potnm(ls, nm) result(out)
      use global, only: prof_nm
      integer :: ls, nm
      real(r8) :: out(nm, ls)

      out = prof_nm

   end function get_potnm


   subroutine test_meshgrid(a, b, la, lb, aa, bb)
      use fft_mod, only: meshgrid
      integer, intent(in) :: la, lb
      real(r8), intent(in) :: a(la), b(lb)
      real(r8), intent(out)  :: aa(lb, la), bb(lb, la)

      call meshgrid(a, b, aa, bb)

   end subroutine test_meshgrid


   subroutine test_gradient(ns, nth, nph, grid, grad)
      use derivatives
      use global, only : len_s, len_ph, len_th
      integer, intent(in) :: ns, nth, nph
      complex(r8), intent(in) :: grid(ns, nth, nph)
      real(r8), intent(out) :: grad(ns, nth, nph, 3)

      grad = gradient(grid)

   end subroutine test_gradient


   subroutine test_cross(a, b, l1, l2, l3, out)
      use helper
      integer :: l1, l2, l3
      real(r8), intent(in) :: a(l1, l2, l3, 3), b(l1, l2, l3, 3)
      real(r8), intent(out) :: out(l1, l2, l3, 3)

      out = cross_product_real(a, b)

   end subroutine test_cross

   subroutine get_j_super(out, len_s, len_th, len_ph, num_modes)
      use global, only : j_super
      integer(i8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(r8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)

      out = j_super

   end subroutine get_j_super

   subroutine get_j_xyz(out, len_s, len_th, len_ph, num_modes)
      use global, only : j_super, es => e_sub_s, eth => e_sub_th, eph => e_sub_ph
      integer(i8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(r8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)
      integer :: j, k

      !$OMP PARALLEL DO PRIVATE(j, k)
      do k=1, num_modes
         do j=1, 3
            out(:,:,:,j,k) = j_super(:,:,:,1,k) * es(:,:,:,j) + &
               j_super(:,:,:,2,k) * eth(:,:,:,j) + &
               j_super(:,:,:,3,k) * eph(:,:,:,j)
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine get_j_xyz

   subroutine get_gradpar_pot_super(out, len_s, len_th, len_ph, num_modes)
      use global, only : gradpar_pot_super
      integer(i8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(r8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)

      out = gradpar_pot_super

   end subroutine get_gradpar_pot_super

   subroutine get_gradpar_pot_xyz(out, len_s, len_th, len_ph, num_modes)
      use global, only : gradpar_pot_super, es => e_sub_s, eth => e_sub_th, eph => e_sub_ph
      integer(i8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(r8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)
      integer :: j, k

      !$OMP PARALLEL DO PRIVATE(j, k)
      do k=1, num_modes
         do j=1, 3
            out(:,:,:,j,k) = gradpar_pot_super(:,:,:,1,k) * es(:,:,:,j) + &
               gradpar_pot_super(:,:,:,2,k) * eth(:,:,:,j) + &
               gradpar_pot_super(:,:,:,3,k) * eph(:,:,:,j)
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine get_gradpar_pot_xyz

end module synthetic_mirnov

