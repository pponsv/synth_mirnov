module potential
   use constants

   implicit none

contains

   ! function gaussian_profile(s, s0, sigma, amp) result(out)
   !    real(r8), intent(in) :: s(:), s0, sigma, amp
   !    real(r8) :: out(size(s))

   !    out = amp * exp(-((s - s0)**2) / sigma)
   ! end function gaussian_profile

   subroutine init_pot(profiles, ms, ns, fs, time)
      use global
      integer(i8), intent(in) :: ms(:), ns(:)
      real(r8), intent(in) :: fs(:), time(:)
      complex(r8), intent(in) :: profiles(:,:)
      integer :: i

      num_modes = size(ms)
      ms_pot = ms
      ns_pot = ns
      ws_pot = 2.*PI*fs
      t = time
      len_t = size(time)
      prof_nm = profiles

      allocate(dpot_dph(len_s, len_th, len_ph, num_modes))
      allocate(dpot_dth(len_s, len_th, len_ph, num_modes))
      allocate(gradpar_pot_super(len_s, len_th, len_ph, 3, num_modes))
      allocate(j_super(len_s, len_th, len_ph, 3, num_modes))

   end subroutine init_pot

   subroutine potential_gradients
      use global
      integer :: i, j, k, l

      !$OMP PARALLEL DO PRIVATE(i, j, k, l)
      do l=1, num_modes
         do k=1, len_ph
            do j=1, len_th
               do i=1, len_s
                  dpot_dth(i, j, k, l) = prof_nm(l,i) * (imag * ms_pot(l)) * exp(imag*(ms_pot(l) * th(j) + ns_pot(l) * ph(k)))
                  dpot_dph(i, j, k, l) = prof_nm(l,i) * (imag * ns_pot(l)) * exp(imag*(ms_pot(l) * th(j) + ns_pot(l) * ph(k)))
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(k)
      do k=1, num_modes
         gradpar_pot_super(:,:,:,1,k) = 0
         gradpar_pot_super(:,:,:,2,k) = inv_mod_b2 * (b_super(:,:,:,2) * dpot_dth(:,:,:,k) + &
            b_super(:,:,:,3) * dpot_dph(:,:,:,k))
         gradpar_pot_super(:,:,:,3,k) = gradpar_pot_super(:,:,:,2,k) * b_super(:,:,:,3)
         gradpar_pot_super(:,:,:,2,k) = gradpar_pot_super(:,:,:,2,k) * b_super(:,:,:,2)
      end do
      !$OMP END PARALLEL DO

   end subroutine potential_gradients

   subroutine potential_curls
      use global
      use helper
      use derivatives
      ! type(vector_grid) :: curl_pot
      complex(r8), allocatable :: curl_pot(:,:,:,:), tmp(:,:,:,:)
      integer :: l

      do l=1, num_modes
         curl_pot = curl(lower_indices(gradpar_pot_super(:,:,:,:,l), g_sub_ij), inv_sqrt_g)
         j_super(:,:,:,:,l) = curl(lower_indices(curl_pot, g_sub_ij), inv_sqrt_g)
      end do

   end subroutine potential_curls

end module potential
