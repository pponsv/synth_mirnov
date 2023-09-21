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
         j_super = curl(lower_indices(curl(gradpar_pot_sub, inv_sqrt_g), g_sub_ij), inv_sqrt_g)

         !Pasar curl a cartesianas (multiplicando por e_sub_i)
         print *, 'TIMESTEP: ', idx_time, '  t=', t(idx_time)

         !$OMP PARALLEL DO PRIVATE(idx_coil, integrand)
         do idx_coil=1, num_coils
            
            integrand = sum_three_vector_grids(&
               product_with_scalar_grid(r_x_e_s(idx_coil),  j_super%u1), &
               product_with_scalar_grid(r_x_e_th(idx_coil), j_super%u2), &
               product_with_scalar_grid(r_x_e_ph(idx_coil), j_super%u3))
            db_coils(idx_coil, 1, idx_time) = sum(integrand%u1)
            db_coils(idx_coil, 2, idx_time) = sum(integrand%u2)
            db_coils(idx_coil, 3, idx_time) = sum(integrand%u3)
         end do
         !$OMP END PARALLEL DO
      end do
   end subroutine main_loop

   subroutine init_coils_main
      use global, only : coil_positions, num_coils, r_coil, xyz_grid, &
         r_3_sqrtg, len_s, len_th, len_ph, sqrt_g, &
         r_x_e_s, r_x_e_th, r_x_e_ph, e_sub_s, e_sub_th, e_sub_ph
      use types
      integer :: idx_coil

      allocate(r_3_sqrtg(num_coils, len_s, len_th, len_ph))
      allocate(r_coil(num_coils), &
         r_x_e_s(num_coils), &
         r_x_e_th(num_coils), &
         r_x_e_ph(num_coils))

      !$OMP PARALLEL DO PRIVATE(idx_coil)
      do idx_coil = 1, num_coils
         r_coil(idx_coil)%u1 = xyz_grid%u1 - coil_positions(1,idx_coil)
         r_coil(idx_coil)%u2 = xyz_grid%u2 - coil_positions(2,idx_coil)
         r_coil(idx_coil)%u3 = xyz_grid%u3 - coil_positions(3,idx_coil)
         r_3_sqrtg(idx_coil, :, :, :) = PI * sqrt_g * sqrt(dot_vector_grid(r_coil(idx_coil), r_coil(idx_coil)))**(-3)
         r_x_e_s(idx_coil)  = product_with_scalar_grid(cross_vector_grid(e_sub_s,  r_coil(idx_coil)), r_3_sqrtg(idx_coil, :, :, :))
         r_x_e_th(idx_coil) = product_with_scalar_grid(cross_vector_grid(e_sub_th, r_coil(idx_coil)), r_3_sqrtg(idx_coil, :, :, :))
         r_x_e_ph(idx_coil) = product_with_scalar_grid(cross_vector_grid(e_sub_ph, r_coil(idx_coil)), r_3_sqrtg(idx_coil, :, :, :))
      end do
      !$OMP END PARALLEL DO
   end subroutine init_coils_main
end module main
