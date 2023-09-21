module main

contains
   subroutine main_loop(db_coils)
      use global
      use potential, only : potential_gradients
      use derivatives, only: curl, gradient
      use types
      real(8), intent(out) :: db_coils(num_coils, 3, len_t)
      integer :: idx_time
      integer :: idx_coil
      integer :: i, j, k
      real(8) :: tmp
      type(vector_grid) :: integrand

      print *, "Main loop", num_coils, size(db_coils, 1), size(db_coils, 2), size(db_coils, 3)

      ! $OMP PARALLEL DO PRIVATE(idx_time, idx_coil, integrand, gradpar_pot_sub, curl_curl_gradpar_pot)
      do idx_time=1, size(t)

         call potential_gradients(t(idx_time))
         gradpar_pot_sub = lower_indices(gradpar_pot_super, g_sub_ij)
         curl_curl_gradpar_pot = curl(lower_indices(curl(gradpar_pot_sub), g_sub_ij))

         !Pasar curl a cartesianas (multiplicando por e_sub_i)
         print *, 'TIMESTEP: ', idx_time, '  t=', t(idx_time)

         ! $OMP PARALLEL DO PRIVATE(idx_coil, integrand)
         do idx_coil=1, num_coils
            ! print *, idx_coil, "debug"
            integrand = cross_vector_grid(curl_curl_gradpar_pot, r_coil(idx_coil))
            ! $OMP WORKSHARE
            db_coils(idx_coil, 1, idx_time) = sum(integrand%u1 * vol_element)
            db_coils(idx_coil, 2, idx_time) = sum(integrand%u2 * vol_element)
            db_coils(idx_coil, 3, idx_time) = sum(integrand%u3 * vol_element)
            ! $OMP END WORKSHARE
         end do
         ! $OMP END PARALLEL DO
      end do
      ! $OMP END PARALLEL DO
   end subroutine main_loop
end module main
