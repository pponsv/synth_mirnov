module potential
   use constants

   implicit none

contains

   subroutine try_dealloc_pot
      use global, only : dpot_dph, dpot_dth, gradpar_pot_super, j_super

      if (allocated(dpot_dph)) deallocate(dpot_dph)
      if (allocated(dpot_dth)) deallocate(dpot_dth)
      if (allocated(gradpar_pot_super)) deallocate(gradpar_pot_super)
      if (allocated(j_super)) deallocate(j_super)

   end subroutine try_dealloc_pot


   subroutine init_pot(profiles, ms, ns, fs, time)
      use global
      integer(i8), intent(in) :: ms(:), ns(:)
      real(r8), intent(in) :: fs(:), time(:)
      complex(r8), intent(in) :: profiles(:,:)

      num_modes = size(ms)
      ms_pot = ms
      ns_pot = ns
      ws_pot = 2.*PI*fs
      t = time
      len_t = size(time)
      prof_nm = profiles

      call try_dealloc_pot

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

      !$OMP PARALLEL DO PRIVATE(i, j, k, l)
      do l=1, num_modes
         do k=1, len_ph
            do j=1, len_th
               do i=1, len_s
                  gradpar_pot_super(i,j,k,1,l) = 0
                  gradpar_pot_super(i,j,k,2,l) = inv_mod_b2(i,j,k) * (b_super(i,j,k,2) * dpot_dth(i,j,k,l) + &
                     b_super(i,j,k,3) * dpot_dph(i,j,k,l))
                  gradpar_pot_super(i,j,k,3,l) = gradpar_pot_super(i,j,k,2,l) * b_super(i,j,k,3)
                  gradpar_pot_super(i,j,k,2,l) = gradpar_pot_super(i,j,k,2,l) * b_super(i,j,k,2)
               end do
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end subroutine potential_gradients


   subroutine potential_curls
      use global
      use helper
      use derivatives
      complex(r8) :: curl_pot(len_s, len_th, len_ph, 3)
      integer :: l

      do l=1, num_modes
         curl_pot = curl(lower_indices(gradpar_pot_super(:,:,:,:,l), g_sub_ij), inv_sqrt_g)
         j_super(:,:,:,:,l) = curl(lower_indices(curl_pot, g_sub_ij), inv_sqrt_g)
      end do

   end subroutine potential_curls

end module potential
