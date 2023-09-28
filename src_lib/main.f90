module main

   use constants

contains

   subroutine main_loop(db_coils)
      use global
      use types
      use potential, only : potential_gradients, potential_curls

      real(r8), intent(out) :: db_coils(3, num_coils, len_t)
      complex(r8) :: db(3, num_modes, num_coils)
      integer :: idx_time, idx_coil, idx_mode
      integer :: i, j, k
      real(r8) :: int_factor
      type(vector_grid) :: us, uth, uph

      print *, "MAIN LOOP"
      print *, "num_coils=", num_coils, "len_t=", len_t
      print *, "GRID SIZE=", len_s, len_th, len_ph
      print *, "DB SIZE=", size(db_coils, 1), size(db_coils, 2), size(db_coils, 3)

      int_factor = delta_s * delta_th * delta_ph
      db_coils = 0.
      db = 0.
      tmp_b = 0.

      ! print *, e_sub_s%u1

      call potential_gradients

      call potential_curls

      do idx_coil=1, num_coils

         print *, 'COIL: ', idx_coil, '  (total ', num_coils, ')'

         call init_coil(idx_coil, us, uth, uph)

         do idx_mode=1, num_modes
            db(:, idx_mode, idx_coil) = integrate_mode(idx_mode, int_factor, us, uth, uph)
         end do

         do idx_time=1, len_t
            do idx_mode=1, num_modes
               db_coils(:, idx_coil, idx_time) = db_coils(:, idx_coil, idx_time) - &
                  db(:, idx_mode, idx_coil)*exp( -imag * ws_pot(idx_mode) * t(idx_time)) / (4*pi)
            end do
         end do

      end do
      ! print *, tmp_b(:,:,1)

   end subroutine main_loop


   function integrate_mode(idx_mode, int_factor, es, eth, eph) result(db)
      use types
      use global, only: len_s, len_th, len_ph, j_super
      type(vector_grid), intent(in) :: es, eth, eph
      real(r8) :: int_factor
      complex(r8) :: dbx, dby, dbz, db(3)
      db = 0.
      dbx = 0.
      dby = 0.
      dbz = 0.

      !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION (+:dbx, dby, dbz)
      do k=1, len_ph
         do j=1, len_th
            do i=1, len_s
               dbx = dbx + (es%u1(i, j, k) * j_super(idx_mode)%u1(i, j, k) + &
                  eth%u1(i, j, k) * j_super(idx_mode)%u2(i, j, k) + &
                  eph%u1(i, j, k) * j_super(idx_mode)%u3(i, j, k)) * int_factor
               dby = dby + (es%u2(i, j, k) * j_super(idx_mode)%u1(i, j, k) + &
                  eth%u2(i, j, k) * j_super(idx_mode)%u2(i, j, k) + &
                  eph%u2(i, j, k) * j_super(idx_mode)%u3(i, j, k)) * int_factor
               dbz = dbz + (es%u3(i, j, k) * j_super(idx_mode)%u1(i, j, k) + &
                  eth%u3(i, j, k) * j_super(idx_mode)%u2(i, j, k) + &
                  eph%u3(i, j, k) * j_super(idx_mode)%u3(i, j, k))* int_factor
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      db = (/ dbx, dby, dbz /)

   end function integrate_mode


   subroutine init_coil(idx_coil, us, uth, uph)
      use types
      ! use global
      use global, only : coil_xyz, xyz_grid, e_sub_s, e_sub_ph, e_sub_th,&
         len_s, len_th, len_ph, sqrt_g
      type(vector_grid), intent(out) :: us, uth, uph
      type(vector_grid) :: r_coil
      real(r8) :: r_3_sqrtg(len_s, len_th, len_ph)
      integer, intent(in) :: idx_coil

      r_coil%u1 = coil_xyz(1,idx_coil) - xyz_grid%u1
      r_coil%u2 = coil_xyz(2,idx_coil) - xyz_grid%u2
      r_coil%u3 = coil_xyz(3,idx_coil) - xyz_grid%u3

      r_3_sqrtg(:, :, :) = sqrt_g * sqrt(dot_vector_grid(r_coil, r_coil))**(-3)

      us = product_with_scalar_grid(cross_vector_grid(e_sub_s,   r_coil), r_3_sqrtg)
      uth = product_with_scalar_grid(cross_vector_grid(e_sub_th, r_coil), r_3_sqrtg)
      uph = product_with_scalar_grid(cross_vector_grid(e_sub_ph, r_coil), r_3_sqrtg)

   end subroutine init_coil

end module main
