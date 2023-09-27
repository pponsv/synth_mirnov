module main

contains
   subroutine main_loop(db_coils)
      use global
      use potential, only : potential_gradients
      use derivatives, only: curl, gradient
      use types
      real(r8), intent(out) :: db_coils(num_coils, 3, len_t)
      integer :: idx_time
      integer :: idx_coil
      integer :: i, j, k
      real(r8) :: dbx, dby, dbz, int_factor
      type(vector_grid) :: tmp

      print *, "MAIN LOOP"
      print *, "num_coils=", num_coils, "len_t=", len_t
      print *, "GRID SIZE=", len_s, len_th, len_ph

      int_factor = delta_s * delta_th * delta_ph
      db_coils = 0.

      do idx_time=1, size(t)

         print *, 'TIMESTEP: ', idx_time, '  t=', t(idx_time)

         call potential_gradients(t(idx_time))

         j_super = curl(lower_indices(gradpar_pot_super, g_sub_ij), inv_sqrt_g)
         j_super = curl(lower_indices(j_super, g_sub_ij), inv_sqrt_g)

         ! $OMP PARALLEL DO PRIVATE(idx_coil, integrand)
         do idx_coil=1, num_coils
            dbx = 0.
            dby = 0.
            dbz = 0.
            !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION (+:dbx, dby, dbz)
            do k=1, len_ph
               do j=1, len_th
                  do i=1, len_s
                     dbx = dbx + (e_sub_s_x_r(idx_coil)%u1(i, j, k) * j_super%u1(i, j, k) + &
                        e_sub_th_x_r(idx_coil)%u1(i, j, k) * j_super%u2(i, j, k) + &
                        e_sub_ph_x_r(idx_coil)%u1(i, j, k) * j_super%u3(i, j, k)) * int_factor
                     dby = dby + (e_sub_s_x_r(idx_coil)%u2(i, j, k) * j_super%u1(i, j, k) + &
                        e_sub_th_x_r(idx_coil)%u2(i, j, k) * j_super%u2(i, j, k) + &
                        e_sub_ph_x_r(idx_coil)%u2(i, j, k) * j_super%u3(i, j, k)) * int_factor
                     dbz = dbz + (e_sub_s_x_r(idx_coil)%u3(i, j, k) * j_super%u1(i, j, k) + &
                        e_sub_th_x_r(idx_coil)%u3(i, j, k) * j_super%u2(i, j, k) + &
                        e_sub_ph_x_r(idx_coil)%u3(i, j, k) * j_super%u3(i, j, k))* int_factor
                  end do
               end do
            end do
            !$OMP END PARALLEL DO
            db_coils(idx_coil, 1, idx_time) = -dbx/(4*pi)
            db_coils(idx_coil, 2, idx_time) = -dby/(4*pi)
            db_coils(idx_coil, 3, idx_time) = -dbz/(4*pi)

            !    integrand = sum_three_vector_grids(&
            !       product_with_scalar_grid(r_x_e_s(idx_coil),  j_super%u1), &
            !       product_with_scalar_grid(r_x_e_th(idx_coil), j_super%u2), &
            !       product_with_scalar_grid(r_x_e_ph(idx_coil), j_super%u3))
            !    db_coils(idx_coil, 1, idx_time) = sum(integrand%u1)
            !    db_coils(idx_coil, 2, idx_time) = sum(integrand%u2)
            !    db_coils(idx_coil, 3, idx_time) = sum(integrand%u3)
         end do
         ! !$OMP END PARALLEL DO
      end do
   end subroutine main_loop

   subroutine init_booz

   end subroutine init_booz

   subroutine init_coils_main
      use global, only : coil_positions, num_coils, r_coil, xyz_grid, &
         r_3_sqrtg, len_s, len_th, len_ph, sqrt_g, &
         e_sub_s_x_r, e_sub_th_x_r, e_sub_ph_x_r, e_sub_s, e_sub_th, e_sub_ph
      use types
      integer :: idx_coil

      allocate(r_3_sqrtg(len_s, len_th, len_ph, num_coils))
      allocate(r_coil(num_coils), &
         e_sub_s_x_r(num_coils), &
         e_sub_th_x_r(num_coils), &
         e_sub_ph_x_r(num_coils))

      !$OMP PARALLEL DO PRIVATE(idx_coil)
      do idx_coil = 1, num_coils
         r_coil(idx_coil)%u1 = coil_positions(1,idx_coil) - xyz_grid%u1
         r_coil(idx_coil)%u2 = coil_positions(2,idx_coil) - xyz_grid%u2
         r_coil(idx_coil)%u3 = coil_positions(3,idx_coil) - xyz_grid%u3
         r_3_sqrtg(:, :, :, idx_coil) = sqrt_g * sqrt(dot_vector_grid(r_coil(idx_coil), r_coil(idx_coil)))**(-3)
         e_sub_s_x_r(idx_coil)  = product_with_scalar_grid(cross_vector_grid(e_sub_s,  r_coil(idx_coil)), &
            r_3_sqrtg(:, :, :, idx_coil))
         e_sub_th_x_r(idx_coil) = product_with_scalar_grid(cross_vector_grid(e_sub_th, r_coil(idx_coil)), &
            r_3_sqrtg(:, :, :, idx_coil))
         e_sub_ph_x_r(idx_coil) = product_with_scalar_grid(cross_vector_grid(e_sub_ph, r_coil(idx_coil)), &
            r_3_sqrtg(:, :, :, idx_coil))
      end do
      !$OMP END PARALLEL DO
   end subroutine init_coils_main
end module main
