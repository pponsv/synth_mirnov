module main

   use constants

   implicit none

contains

   subroutine main_loop(db_coils)
      use global
      use potential, only : potential_gradients, potential_curls

      real(r8), intent(out) :: db_coils(3, num_coils, len_t)
      complex(r8) :: db(3, num_modes, num_coils)
      integer :: idx_time, idx_coil, idx_mode
      real(r8) :: int_factor
      real(r8), dimension(len_s, len_th, len_ph, 3)  :: us, uth, uph
      integer :: i, j, k

      write (*, '(/, A, /)') "MAIN LOOP"
      write (*, '(2(A12, I5, 3X, /))') "NUM COILS = ", num_coils, "len_t =     ", len_t
      write (*, '(A12, 3(I5, 3X))') "GRID SIZE = ", len_s, len_th, len_ph
      write (*, '(A12, 3(I5, 3X))') "DB SIZE =   ", size(db_coils, 1), size(db_coils, 2), size(db_coils, 3)

      int_factor = delta_s * delta_th * delta_ph
      db_coils = 0.
      db = 0.

      call potential_gradients

      call potential_curls

      do idx_coil=1, num_coils

         write (*, '(A6, I3, A3, I0)') 'COIL: ', idx_coil, ' / ', num_coils

         call init_coil(idx_coil, us, uth, uph)

         do idx_mode=1, num_modes
            db(:, idx_mode, idx_coil) = integrate_mode(idx_mode, int_factor, us, uth, uph)
         end do

         do idx_time=1, len_t
            do idx_mode=1, num_modes
               db_coils(:, idx_coil, idx_time) = db_coils(:, idx_coil, idx_time) + &
                  db(:, idx_mode, idx_coil)*exp( -imag * ws_pot(idx_mode) * t(idx_time))
            end do
         end do

      end do
      ! print *, tmp_b(:,:,1)

   end subroutine main_loop


   function integrate_mode(idx_mode, int_factor, es, eth, eph) result(db)
      ! use types
      use global, only: len_s, len_th, len_ph, j_super
      integer, intent(in) :: idx_mode
      real(r8), intent(in), dimension(len_s, len_th, len_ph, 3) :: es, eth, eph
      real(r8) :: int_factor
      complex(r8) :: dbx, dby, dbz, db(3)
      integer :: i, j, k
      db = 0.
      dbx = 0.
      dby = 0.
      dbz = 0.

      !$OMP PARALLEL DO PRIVATE(i, j, k) REDUCTION (+:dbx, dby, dbz)
      do k=1, len_ph
         do j=1, len_th
            do i=1, len_s
               dbx = dbx + (es(i, j, k, 1) * j_super(i, j, k, 1, idx_mode) + &
                  eth(i, j, k, 1) * j_super(i, j, k, 2, idx_mode) + &
                  eph(i, j, k, 1) * j_super(i, j, k, 3, idx_mode))
               dby = dby + (es(i, j, k, 2) * j_super(i, j, k, 1, idx_mode) + &
                  eth(i, j, k, 2) * j_super(i, j, k, 2, idx_mode) + &
                  eph(i, j, k, 2) * j_super(i, j, k, 3, idx_mode))
               dbz = dbz + (es(i, j, k, 3) * j_super(i, j, k, 1, idx_mode) + &
                  eth(i, j, k, 3) * j_super(i, j, k, 2, idx_mode) + &
                  eph(i, j, k, 3) * j_super(i, j, k, 3, idx_mode))
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      db = - (/ dbx, dby, dbz /) * int_factor / (4*pi)

   end function integrate_mode


   subroutine init_coil(idx_coil, us, uth, uph)
      use helper
      use global, only : coil_xyz, xyz_grid, e_sub_s, e_sub_ph, e_sub_th,&
         len_s, len_th, len_ph, sqrt_g
      real(r8), intent(out), dimension(len_s, len_th, len_ph, 3) :: us, uth, uph
      real(r8), dimension(len_s, len_th, len_ph, 3) :: r_coil
      real(r8) :: r_3_sqrtg(len_s, len_th, len_ph)
      integer, intent(in) :: idx_coil

      r_coil(:,:,:,1) = coil_xyz(1,idx_coil) - xyz_grid(:,:,:,1)
      r_coil(:,:,:,2) = coil_xyz(2,idx_coil) - xyz_grid(:,:,:,2)
      r_coil(:,:,:,3) = coil_xyz(3,idx_coil) - xyz_grid(:,:,:,3)

      r_3_sqrtg = sqrt_g * sqrt(dot_product_real(r_coil, r_coil))**(-3)

      us  = scalar_product_real(cross_product_real(e_sub_s,  r_coil), r_3_sqrtg)
      uth = scalar_product_real(cross_product_real(e_sub_th, r_coil), r_3_sqrtg)
      uph = scalar_product_real(cross_product_real(e_sub_ph, r_coil), r_3_sqrtg)

   end subroutine init_coil

end module main
