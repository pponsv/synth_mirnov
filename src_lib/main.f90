module main

   use constants

   implicit none

contains

   subroutine main_loop(db_coils)
      use global
      use potential, only : potential_gradients, potential_curls

      complex(r8), intent(out) :: db_coils(3, num_coils, len_t)

      write (*, '(/, A)', advance="no") "MAIN LOOP: gradients... "
      call potential_gradients

      write (*, '(/, A)', advance="no") "curls... "
      call potential_curls

      write (*, '(/, A)', advance="no") "loop ..."
      call loop_over_coils(db_coils)

   end subroutine main_loop


   subroutine loop_over_coils(db_coils)
      use global
      integer :: idx_time, idx_coil, idx_mode
      complex(r8), intent(out) :: db_coils(3, num_coils, len_t)
      complex(r8) :: db(3, num_modes, num_coils)

      db_coils = 0.
      db = 0.

      !!$OMP PARALLEL DO PRIVATE(idx_coil, idx_mode)
      do idx_coil=1, num_coils

         !  Integral for each mode
         do idx_mode=1, num_modes
            db(:, idx_mode, idx_coil) = integrate_mode(idx_mode, idx_coil, int_factor)
         end do
      end do
      !!$OMP END PARALLEL DO

      do idx_coil=1, num_coils
         !  Time evolution
         do idx_time=1, len_t
            do idx_mode=1, num_modes
               db_coils(:, idx_coil, idx_time) = db_coils(:, idx_coil, idx_time) + &
                  db(:, idx_mode, idx_coil)*exp( -IMAG * ws_pot(idx_mode) * t(idx_time))
            end do
         end do

      end do
   end subroutine loop_over_coils


   function integrate_mode(idx_mode, idx_coil, int_factor) result(db_int)
      use global, only : len_s, len_th, len_ph, j_super, us, uth, uph, db_all
      use helper, only : scalar_product_real_cplx
      integer, intent(in) :: idx_mode, idx_coil
      real(r8), intent(in) :: int_factor
      complex(r8) :: dbx, dby, dbz, db_int(3)

      integer :: i, j, k
      db_int = 0.
      dbx = 0.
      dby = 0.
      dbz = 0.

      db_all(:,:,:,:,idx_mode) = (&
         scalar_product_real_cplx(us(:,:,:,:,idx_coil), j_super(:,:,:,1,idx_mode)) + &
         scalar_product_real_cplx(uth(:,:,:,:,idx_coil), j_super(:,:,:,2,idx_mode)) + &
         scalar_product_real_cplx(uph(:,:,:,:,idx_coil), j_super(:,:,:,3,idx_mode)))

      do k=1, len_ph
         do j=1, len_th
            do i=1, len_s
               dbx = dbx + db_all(i, j, k, 1, idx_mode)
               dby = dby + db_all(i, j, k, 2, idx_mode)
               dbz = dbz + db_all(i, j, k, 3, idx_mode)
            end do
         end do
      end do

      db_int = - (/ dbx, dby, dbz /) * int_factor / (4*pi)

   end function integrate_mode


end module main
