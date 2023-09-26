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

      delta_s = s(2) - s(1)
      delta_th = th(2) - th(1)
      delta_ph = ph(2) - ph(1)

      th_freqs = 2*PI*fftfreqs(len_th, delta_th)
      ph_freqs = 2*PI*fftfreqs(len_ph, delta_ph)
      allocate(ftheta(len_th, len_ph), fphi(len_th, len_ph))
      call meshgrid(ph_freqs, th_freqs, fphi, ftheta)

      call plan_ffts()

      iota = iota_b
      mod_b = b_mod_b
      inv_mod_b2 = 1. / b_mod_b**2
      sqrt_g = sqrt_g_b
      inv_sqrt_g = 1. / sqrt_g_b

      call b_super%alloc(len_s, len_th, len_ph)
      ! allocate(B_super%u1, B_super%u2, B_super%u3, mold=b_mod_b)
      B_super%u1 = 0
      B_super%u3 = -1/(2*pi*sqrt_g)*phi_b_g ! 2*pi sobra??
      ! B_super%u3 = -1/(sqrt_g)*phi_b_g
      do i=1, len_s
         B_super%u2(i, :, :) = B_super%u3(i, :, :) * iota(i)
      end do

      xyz_grid%u1 = x
      xyz_grid%u2 = y
      xyz_grid%u3 = z

      call gradpar_pot_super%alloc(len_s, len_ph, len_th)

   end subroutine init_booz

   subroutine init_basis(e_sub_s_b, e_sub_th_b, e_sub_ph_b, len_s_b, len_th_b, len_ph_b)
      use global, only : e_sub_s, e_sub_th, e_sub_ph, g_sub_ij, b_super
      integer(8), intent(in) :: len_s_b, len_ph_b, len_th_b
      real(8), intent(in), dimension(3, len_s_b, len_th_b, len_ph_b) :: e_sub_s_b, e_sub_th_b, e_sub_ph_b
      integer :: i, j, k

      call e_sub_s%init_grid(e_sub_s_b)
      call e_sub_th%init_grid(e_sub_th_b)
      call e_sub_ph%init_grid(e_sub_ph_b)

      call g_sub_ij%init(e_sub_s, e_sub_th, e_sub_ph)

   end subroutine init_basis

   subroutine init_pot(profiles, ms, ns, fs, time, num_modes, len_time, len_s)
      use potential, only : initialize_potential
      integer(8), intent(in) :: num_modes, len_time, len_s
      integer(8), intent(in) :: ms(num_modes), ns(num_modes)
      real(8), intent(in) :: fs(num_modes), time(len_time)
      complex(8), intent(in) :: profiles(num_modes, len_s)

      call initialize_potential(profiles=profiles, ms=ms, ns=ns, fs=fs, time=time)

   end subroutine init_pot

   subroutine init_coils(xyzs, ncoils)
      use main, only: init_coils_main
      use global, only : coil_positions, num_coils
      integer(8), intent(in) :: ncoils
      real(8), intent(in) :: xyzs(3, ncoils)

      coil_positions = xyzs
      num_coils = ncoils
      call init_coils_main

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

   subroutine test_magnetic_field(i, j, k)
      use global
      integer, intent(in) :: i, j, k
      print *,  b_super%u2(i, j, k) *e_sub_th%u1(i, j, k) + b_super%u3(i, j, k) * e_sub_ph%u1(i, j, k)
      print *,  b_super%u2(i, j, k) *e_sub_th%u2(i, j, k) + b_super%u3(i, j, k) * e_sub_ph%u2(i, j, k)
      print *,  b_super%u2(i, j, k) *e_sub_th%u3(i, j, k) + b_super%u3(i, j, k) * e_sub_ph%u3(i, j, k)
   end subroutine test_magnetic_field

   function get_potnm(ls, nm) result(out)
      use global, only: potnm
      integer :: ls, nm
      real(8) :: out(nm, ls)

      out = potnm
   end function get_potnm

   subroutine test_meshgrid(a, b, la, lb, aa, bb)
      use fft_mod, only: meshgrid
      use global, only: ftheta, fphi
      integer, intent(in) :: la, lb
      real(8), intent(in) :: a(la), b(lb)
      real(8), intent(out)  :: aa(lb, la), bb(lb, la)

      call meshgrid(a, b, aa, bb)
      aa = ftheta
      bb = fphi

      ! print *, 'AA'
      ! print *, aa
      ! print *, 'BB'
      ! print *, bb(1,:)
   end subroutine test_meshgrid

   subroutine test_gradient(ns, nth, nph, grid, grad)
      use Derivatives
      use types
      use global, only : len_s, len_ph, len_th
      integer, intent(in) :: ns, nth, nph
      real(8), intent(in) :: grid(ns, nth, nph)
      type(vector_grid) :: tgrad
      real(8), intent(out) :: grad(3, ns, nth, nph)

      len_s =  ns
      len_ph = nph
      len_th = nth

      tgrad = gradient(grid)

      grad(1,:,:,:) = tgrad%u1
      grad(2,:,:,:) = tgrad%u2
      grad(3,:,:,:) = tgrad%u3

   end subroutine test_gradient


end module synthetic_mirnov
