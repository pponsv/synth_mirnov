module helper

   use constants

   implicit none

contains

   function cross_product_cplx(a, b) result(out)
      complex(r8), dimension(:,:,:,:), intent(in) :: a, b
      complex(r8) :: out(size(a, 1), size(a, 2), size(a, 3), 3)
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(a, 3)
         do j=1, size(a, 2)
            do i=1, size(a, 1)
               out(i, j, k, 1) = a(i, j, k, 2) * b(i, j, k, 3) - b(i, j, k, 2) * a(i, j, k, 3)
               out(i, j, k, 2) = a(i, j, k, 3) * b(i, j, k, 1) - b(i, j, k, 3) * a(i, j, k, 1)
               out(i, j, k, 3) = a(i, j, k, 1) * b(i, j, k, 2) - b(i, j, k, 1) * a(i, j, k, 2)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function cross_product_cplx


   function cross_product_real(a, b) result(out)
      real(r8), dimension(:,:,:,:), intent(in) :: a, b
      real(r8) :: out(size(a, 1), size(a, 2), size(a, 3), 3)
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(a, 3)
         do j=1, size(a, 2)
            do i=1, size(a, 1)
               out(i, j, k, 1) = a(i, j, k, 2) * b(i, j, k, 3) - b(i, j, k, 2) * a(i, j, k, 3)
               out(i, j, k, 2) = a(i, j, k, 3) * b(i, j, k, 1) - b(i, j, k, 3) * a(i, j, k, 1)
               out(i, j, k, 3) = a(i, j, k, 1) * b(i, j, k, 2) - b(i, j, k, 1) * a(i, j, k, 2)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function cross_product_real


   function dot_product_cplx(a, b) result(out)
      complex(r8), dimension(:,:,:,:), intent(in) :: a, b
      complex(r8) :: out(size(a, 1), size(a, 2), size(a, 3))
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(a, 3)
         do j=1, size(a, 2)
            do i=1, size(a, 1)
               out(i, j, k) = sqrt(a(i, j, k, 1) * b(i, j, k, 1) + a(i, j, k, 2) * b(i, j, k, 2) + &
                  a(i, j, k, 3) * b(i, j, k, 3))
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function dot_product_cplx


   function dot_product_real(a, b) result(out)
      real(r8), dimension(:,:,:,:), intent(in) :: a, b
      real(r8) :: out(size(a, 1), size(a, 2), size(a, 3))
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(a, 3)
         do j=1, size(a, 2)
            do i=1, size(a, 1)
               out(i, j, k) = sum(a(i, j, k, :) * b(i, j, k, :))
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function dot_product_real


   function scalar_product_cplx(vec, sca) result(out)
      complex(r8), intent(in) :: vec(:,:,:,:) , sca(:,:,:)
      complex(r8) :: out(size(vec, 1), size(vec, 2), size(vec, 3), 3)
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(vec, 3)
         do j=1, size(vec, 2)
            do i=1, size(vec, 1)
               out(i, j, k, :) = vec(i, j, k, :) * sca(i, j, k)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function scalar_product_cplx


   function scalar_product_real_cplx(vec, sca) result(out)
      real(r8), intent(in) :: vec(:,:,:,:)
      complex(r8), intent(in) :: sca(:,:,:)
      complex(r8) :: out(size(vec, 1), size(vec, 2), size(vec, 3), 3)
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(vec, 3)
         do j=1, size(vec, 2)
            do i=1, size(vec, 1)
               out(i, j, k, :) = vec(i, j, k, :) * sca(i, j, k)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function scalar_product_real_cplx

   function scalar_product_real(vec, sca) result(out)
      real(r8), intent(in) :: vec(:,:,:,:)
      real(r8), intent(in) :: sca(:,:,:)
      real(r8) :: out(size(vec, 1), size(vec, 2), size(vec, 3), 3)
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(vec, 3)
         do j=1, size(vec, 2)
            do i=1, size(vec, 1)
               out(i, j, k, :) = vec(i, j, k, :) * sca(i, j, k)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function scalar_product_real


   function lower_indices(vec, gij) result(out)
      complex(r8), intent(in) :: vec(:,:,:,:)
      real(r8), intent(in) :: gij(:,:,:,:,:)
      complex(r8) :: out(size(vec, 1), size(vec, 2), size(vec, 3), 3)
      integer :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k=1, size(vec, 3)
         do j=1, size(vec, 2)
            do i=1, size(vec, 1)
               out(i, j, k, 1) = gij(i, j, k, 1, 1) * vec(i, j, k, 1) + &
                  gij(i, j, k, 1, 2) * vec(i, j, k, 2) + gij(i, j, k, 1, 3) * vec(i, j, k, 3)

               out(i, j, k, 2) = gij(i, j, k, 2, 1) * vec(i, j, k, 1) + &
                  gij(i, j, k, 2, 2) * vec(i, j, k, 2) + gij(i, j, k, 2, 3) * vec(i, j, k, 3)

               out(i, j, k, 3) = gij(i, j, k, 3, 1) * vec(i, j, k, 1) + &
                  gij(i, j, k, 3, 2) * vec(i, j, k, 2) + gij(i, j, k, 3, 3) * vec(i, j, k, 3)
            end do
         end do
      end do
      !$OMP END PARALLEL DO

   end function lower_indices


   function metric_tensor(es, eth, eph) result(out)
      real(r8), intent(in) :: es(:,:,:,:), eth(:,:,:,:), eph(:,:,:,:)
      real(r8) :: out(size(es, 1), size(es, 2), size(es, 3), 3, 3)

      out(:,:,:,1,1) = dot_product_real(es, es)
      out(:,:,:,1,2) = dot_product_real(es, eth)
      out(:,:,:,1,3) = dot_product_real(es, eph)

      out(:,:,:,2,1) = dot_product_real(eth, es)
      out(:,:,:,2,2) = dot_product_real(eth, eth)
      out(:,:,:,2,3) = dot_product_real(eth, eph)

      out(:,:,:,3,1) = dot_product_real(eph, es)
      out(:,:,:,3,2) = dot_product_real(eph, eth)
      out(:,:,:,3,3) = dot_product_real(eph, eph)

   end function metric_tensor


   subroutine checknans(in, len_in)
      real(r8) :: in(*)
      integer :: i, len_in

      do i=1, len_in
         if (in(i) /= in(i)) then
            print *, 'NaN', in(i), i
            stop
         end if
      end do
      print *, 'No NaNS'
   end subroutine checknans

end module helper
