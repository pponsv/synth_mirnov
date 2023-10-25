module synthetic_mirnov

   implicit none

contains

   subroutine init_booz(s_b, th_b, ph_b, b_mod_b, sqrt_g_b, phi_b_g, iota_b, x, y, z, &
      len_s_b, len_th_b, len_ph_b)
      use initialization, only : initialize_boozer
      integer(8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(8), intent(in) :: s_b(len_s_b), th_b(len_th_b), ph_b(len_ph_b), phi_b_g, iota_b(len_s_b)
      real(8), intent(in), dimension(len_s_b, len_th_b, len_ph_b) :: b_mod_b, sqrt_g_b, x, y, z

      call initialize_boozer(s_b, th_b, ph_b, phi_b_g, iota_b, b_mod_b, sqrt_g_b, x, y, z)

   end subroutine init_booz


   subroutine init_basis(e_sub_s_b, e_sub_th_b, e_sub_ph_b, len_s_b, len_th_b, len_ph_b)
      use initialization, only : initialize_basis
      integer(8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(8), intent(in), dimension(len_s_b, len_th_b, len_ph_b, 3) :: e_sub_s_b, e_sub_th_b, e_sub_ph_b

      call initialize_basis(e_sub_s_b, e_sub_th_b, e_sub_ph_b)

   end subroutine init_basis


   subroutine init_pot(profiles, ms, ns, fs, time, num_modes, len_time, len_s)
      use initialization, only : initialize_potentials
      integer(8), intent(in) :: num_modes, len_time, len_s
      integer(8), intent(in) :: ms(num_modes), ns(num_modes)
      real(8), intent(in) :: fs(num_modes), time(len_time)
      complex(8), intent(in) :: profiles(num_modes, len_s)

      call initialize_potentials(profiles=profiles, ms=ms, ns=ns, fs=fs, time=time)

   end subroutine init_pot


   subroutine init_coils(xyzs, ncoils)
      ! use global, only : coil_xyz, num_coils, len_s, len_th, len_ph, r_coil
      use initialization, only : initialize_coils
      integer(8), intent(in) :: ncoils
      real(8), intent(in) :: xyzs(ncoils, 3)

      call initialize_coils(xyzs, ncoils)

   end subroutine init_coils


   subroutine run(num_coils_tmp, len_t_tmp, db_coils)
      use main
      integer(8), intent(in) :: num_coils_tmp, len_t_tmp
      complex(8), intent(out) :: db_coils(3, num_coils_tmp, len_t_tmp)

      call main_loop(db_coils)

   end subroutine run

   !  TEST SUBROUTINES ETC

   subroutine run_coil_loop(num_coils_tmp, len_t_tmp, db_coils)
      use main
      integer(8), intent(in) :: num_coils_tmp, len_t_tmp
      complex(8), intent(out) :: db_coils(3, num_coils_tmp, len_t_tmp)

      call loop_over_coils(db_coils)

   end subroutine run_coil_loop


   subroutine set_j_super(j_super_b, len_s, len_th, len_ph, num_modes)
      use global, only : j_super
      integer(8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(8), intent(in) :: j_super_b(len_s, len_th, len_ph, 3, num_modes)

      j_super = j_super_b

   end subroutine set_j_super


   function test_fft_main(in) result(out)
      use fft_mod
      complex(8) :: in(:,:), nin(size(in,1), size(in, 2))
      complex(8) :: out(size(in,1), size(in, 2)), cout(size(in,1), size(in, 2))

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
      complex(8) :: out(nm, ls)

      out = prof_nm

   end function get_potnm


   subroutine test_meshgrid(a, b, la, lb, aa, bb)
      use fft_mod, only: meshgrid
      integer, intent(in) :: la, lb
      real(8), intent(in) :: a(la), b(lb)
      real(8), intent(out)  :: aa(lb, la), bb(lb, la)

      call meshgrid(a, b, aa, bb)

   end subroutine test_meshgrid


   subroutine test_gradient(ns, nth, nph, grid, grad)
      use derivatives
      integer, intent(in) :: ns, nth, nph
      complex(8), intent(in) :: grid(ns, nth, nph)
      complex(8), intent(out) :: grad(ns, nth, nph, 3)

      grad = gradient(grid)

   end subroutine test_gradient


   subroutine test_cross(a, b, l1, l2, l3, out)
      use helper
      integer :: l1, l2, l3
      real(8), intent(in) :: a(l1, l2, l3, 3), b(l1, l2, l3, 3)
      real(8), intent(out) :: out(l1, l2, l3, 3)

      out = cross_product_real(a, b)

   end subroutine test_cross


   subroutine get_j_super(out, len_s, len_th, len_ph, num_modes)
      use global, only : j_super, ws_pot
      use constants
      integer(8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)
      complex(8) :: fac(size(ws_pot))
      integer :: k

      fac = - imag / (mu0 * 1000 * ws_pot)
      !$OMP PARALLEL DO PRIVATE(k)
      do k=1, num_modes
         out(:,:,:,:,k) = fac(k) * j_super(:,:,:,:,k)
      end do
      !$OMP END PARALLEL DO

   end subroutine get_j_super


   subroutine get_j_xyz(out, len_s, len_th, len_ph, num_modes)
      use global, only : j_super, es => e_sub_s, eth => e_sub_th, eph => e_sub_ph, ws_pot
      use constants
      integer(8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)
      complex(8) :: fac(size(ws_pot))
      integer :: j, k

      fac = - imag / (mu0 * 1000 * ws_pot)
      !$OMP PARALLEL DO PRIVATE(j, k)
      do k=1, num_modes
         do j=1, 3
            out(:,:,:,j,k) = fac(k) * (j_super(:,:,:,1,k) * es(:,:,:,j) + &
               j_super(:,:,:,2,k) * eth(:,:,:,j) + &
               j_super(:,:,:,3,k) * eph(:,:,:,j))
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine get_j_xyz


   subroutine get_gradpar_pot_super(out, len_s, len_th, len_ph, num_modes)
      use global, only : gradpar_pot_super
      integer(8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)

      out = gradpar_pot_super

   end subroutine get_gradpar_pot_super


   subroutine get_gradpar_pot_xyz(out, len_s, len_th, len_ph, num_modes)
      use global, only : gradpar_pot_super, es => e_sub_s, eth => e_sub_th, eph => e_sub_ph
      integer(8), intent(in) :: len_s, len_th, len_ph, num_modes
      complex(8), intent(out) :: out(len_s, len_th, len_ph, 3, num_modes)
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


   subroutine get_r_3_sqrtg(out, len_s, len_th, len_ph)
      use global, only : r_3_sqrtg
      integer(8), intent(in) :: len_s, len_th, len_ph
      real(8), intent(out) :: out(len_s, len_th, len_ph)

      out = r_3_sqrtg

   end subroutine get_r_3_sqrtg


   subroutine test_findiff(y, dx, lx, order, dy)
      use derivatives, only : fd1 => finite_differences_1, fd2 => finite_differences_2, &
         fd3 => finite_differences_3, fd4 => finite_differences_4
      integer, intent(in) :: lx, order
      real(8), intent(in) :: dx
      complex(8), intent(in) :: y(lx)
      complex(8), intent(out) :: dy(lx)

      if (order .eq. 1) dy = fd1(y, dx)
      if (order .eq. 2) dy = fd2(y, dx)
      if (order .eq. 3) dy = fd3(y, dx)
      if (order .eq. 4) dy = fd4(y, dx)

   end subroutine test_findiff

end module synthetic_mirnov

