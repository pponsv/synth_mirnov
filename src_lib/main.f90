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

      print *, "MAIN LOOP"
      print *, "num_coils=", num_coils, "len_t=", len_t

      do idx_time=1, size(t)

         call potential_gradients(t(idx_time))
         gradpar_pot_sub = lower_indices(gradpar_pot_super, g_sub_ij)
         j_super = curl(lower_indices(curl(gradpar_pot_sub), g_sub_ij))

         !Pasar curl a cartesianas (multiplicando por e_sub_i)
         print *, 'TIMESTEP: ', idx_time, '  t=', t(idx_time)

         do idx_coil=1, num_coils
            ! print *, idx_coil, "debug"
            integrand = cross_vector_grid(j_super, r_coil(idx_coil))
            db_coils(idx_coil, 1, idx_time) = sum(integrand%u1 * vol_element)
            db_coils(idx_coil, 2, idx_time) = sum(integrand%u2 * vol_element)
            db_coils(idx_coil, 3, idx_time) = sum(integrand%u3 * vol_element)
         end do
      end do
   end subroutine main_loop

   subroutine init_coils_main
      use global, only : coil_positions, num_coils, r_coil, xyz_grid, abs_r_coil, &
         inv_abs_r_coil_squared, len_s, len_th, len_ph
      use types, only : dot_vector_grid
      integer :: idx_coil

      allocate(inv_abs_r_coil_squared(num_coils, len_s, len_th, len_ph))
      allocate(r_coil(num_coils))

      ! $OMP PARALLEL DO PRIVATE(idx_coil)
      do idx_coil = 1, num_coils
         r_coil(idx_coil)%u1 = xyz_grid%u1 - coil_positions(1,idx_coil)
         r_coil(idx_coil)%u2 = xyz_grid%u2 - coil_positions(2,idx_coil)
         r_coil(idx_coil)%u3 = xyz_grid%u3 - coil_positions(3,idx_coil)
         inv_abs_r_coil_squared(idx_coil, :, :, :) = dot_vector_grid(r_coil(idx_coil), r_coil(idx_coil))
      end do
      ! $OMP END PARALLEL DO
      abs_r_coil = sqrt(inv_abs_r_coil_squared)
      inv_abs_r_coil_squared = 1./inv_abs_r_coil_squared
      print *, (sum(abs_r_coil))
   end subroutine init_coils_main
end module main
