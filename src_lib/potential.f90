module potential
   use constants
   use helper
   implicit none

contains

   ! function gaussian_profile(s, s0, sigma, amp) result(out)
   !    real(8), intent(in) :: s(:), s0, sigma, amp
   !    real(8) :: out(size(s))

   !    out = amp * exp(-((s - s0)**2) / sigma)
   ! end function gaussian_profile

   subroutine initialize_potential(profiles, ms, ns, fs, time)
      use global
      integer(8), intent(in) :: ms(:), ns(:)
      real(8), intent(in) :: fs(:), time(:)
      complex(8) :: profiles(:,:)
      integer :: i

      num_modes = size(ms)
      ms_pot = ms
      ns_pot = ns
      ws_pot = 2.*PI*fs
      t = time
      len_t = size(time)
      potnm = profiles

      exp_mth_pot = exp(complex(0, 1) * outer(real(ms_pot, 8), th))
      exp_nph_pot = exp(complex(0, 1) * outer(real(ns_pot, 8), ph))

      allocate(dpot_dph(len_s, len_th, len_ph), &
         dpot_dth(len_s, len_th, len_ph))

   end subroutine initialize_potential

   subroutine potential_gradients(time)
      use global
      real(8), intent(in) :: time
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, len_ph
         do j=1, len_th
            do i=1, len_s
               dpot_dth(i, j, k) = sum(potnm(:,i) * (complex(0,1) * ms_pot) * &
                  exp_mth_pot(:,j) * exp_nph_pot(:,k) * exp(complex(0, -1) * ws_pot))
               dpot_dph(i, j, k) = sum(potnm(:,i) * (complex(0,1) * ns_pot) * &
                  exp_mth_pot(:,j) * exp_nph_pot(:,k) * exp(complex(0, -1) * ws_pot))
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      gradpar_pot_super%u1 = 0
      gradpar_pot_super%u2 = inv_mod_b2 * (b_super%u2 * dpot_dth + b_super%u3 * dpot_dph)
      gradpar_pot_super%u3 = gradpar_pot_super%u2 * b_super%u3
      gradpar_pot_super%u2 = gradpar_pot_super%u2 * b_super%u2 !  Lo hago así para hacer menos operaciones

   end subroutine potential_gradients


end module potential
