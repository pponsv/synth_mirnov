module potential
   use constants

   implicit none

contains

   subroutine potential_gradients
      use global
      integer :: i, j, k, l

      !$OMP PARALLEL DO PRIVATE(i, j, k, l)
      do l=1, num_modes
         do k=1, len_ph
            do j=1, len_th
               do i=1, len_s
                  dpot_dth(i, j, k, l) = prof_nm(l,i) * (IMAG * ms_pot(l)) * exp(IMAG*(ms_pot(l) * th(j) + ns_pot(l) * ph(k)))
                  dpot_dph(i, j, k, l) = prof_nm(l,i) * (IMAG * ns_pot(l)) * exp(IMAG*(ms_pot(l) * th(j) + ns_pot(l) * ph(k)))
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
