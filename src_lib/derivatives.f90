module Derivatives

   use fft_mod
   ! use helper
   use types

contains


   function gradient(in) result(out)
      use global, only : s, len_s, len_ph, len_th, delta_s, ftheta, fphi
      real(8), intent(in) :: in(len_s, len_th, len_ph)
      type(vector_grid) :: out
      integer :: i, j, k

      call out%alloc(len_s, len_th, len_ph)

      !$OMP PARALLEL DO PRIVATE(k, j)
      do k=1, len_ph
         do j=1, len_th
            out%u1(:, j, k) = finite_differences(in(:, j, k), delta_s)
         end do
      end do
      !$OMP END PARALLEL DO

      !$OMP PARALLEL DO PRIVATE(i)
      do i=1, len_s
         out%u2(i, :, :) = ifft_2d(complex(0,1) * ftheta * fft_2d(in(i, :, :)))
         out%u3(i, :, :) = ifft_2d(complex(0,1) * fphi * fft_2d(in(i, :, :)))
      end do
      !$OMP END PARALLEL DO
   end function gradient

   function curl(in, inv_sqrtg) result(out)
      real(8), intent(in) :: inv_sqrtg(:,:,:)
      type(vector_grid), intent(in):: in
      type(vector_grid) :: out, grad_1, grad_2, grad_3

      ! call out%alloc(size(inv_sqrtg, 1, kind=8), size(inv_sqrtg, 2, kind=8), size(inv_sqrtg, 3, kind=8))

      grad_1 = gradient(in%u1)
      grad_2 = gradient(in%u2)
      grad_3 = gradient(in%u3)

      ! $OMP WORKSHARE
      out%u1 = (grad_3%u2 - grad_2%u3) * inv_sqrtg
      out%u2 = (grad_1%u3 - grad_3%u1) * inv_sqrtg
      out%u3 = (grad_2%u1 - grad_1%u2) * inv_sqrtg
      ! $OMP END WORKSHARE

   end function curl

   pure function finite_differences(y, dx) result(dy)
      real(8), intent(in) :: y(:), dx
      real(8) :: dy(size(y))
      integer :: i, len_y

      len_y = size(y)

      dy(1) = (-3 * y(1) + 4 * y(2) - 1 * y(3)) / (2 * dx)
      dy(2) = (-2 * y(1) - 3 * y(2) + 6 * y(3) - 1 * y(4)) / (6 * dx)
      do i=3, size(y)-2
         dy(i) = (1 * y(i - 2) - 8 * y(i - 1) + 8 * y(i + 1) - 1 * y(i + 2)) / (12 * dx)
      end do
      dy(size(y)-1) = (1 * y(len_y-3) - 6 * y(len_y-2) + 3 * y(len_y-1) + 2 * y(len_y)) / (6 * dx)
      dy(size(y)) = (1 * y(len_y-2) - 4 * y(len_y-1) + 3 * y(len_y)) / (2 * dx)

   end function finite_differences
end module Derivatives
