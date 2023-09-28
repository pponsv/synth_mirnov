module derivatives

   use constants
   use fft_mod
   ! use helper
   use types

contains


   function gradient(in) result(out)
      use global, only : s, len_s, len_ph, len_th, delta_s, &
         delta_th, delta_ph, fth, fph
      complex(r8), intent(in) :: in(len_s, len_th, len_ph)
      type(vector_grid) :: out
      integer :: i, j, k

      call out%alloc(len_s, len_th, len_ph)

      out%u1 = partial_rad(in, delta_s)
      out%u2 = partial_fft(in, fth)
      out%u3 = partial_fft(in, fph)

   end function gradient

   function partial_rad(in, ds) result(out)
      complex(8), intent(in) :: in(:,:,:)
      real(8) :: ds
      complex(8) :: out(size(in, 1), size(in, 2), size(in, 3))
      integer :: j, k

      !$OMP PARALLEL DO PRIVATE(k, j)
      do k=1, size(in, 3)
         do j=1, size(in, 2)
            out(:, j, k) = finite_differences(in(:, j, k), ds)
         end do
      end do
      !$OMP END PARALLEL DO

   end function partial_rad

   function partial_fft(in, coef_grid) result(out)
      complex(8), intent(in) :: in(:,:,:)
      real(r8) :: coef_grid(:,:)
      complex(8) :: out(size(in, 1), size(in, 2), size(in, 3))

      !$OMP PARALLEL DO PRIVATE(i)
      do i=1, size(in, 1)
         out(i, :, :) = ifft_2d(complex(0,1) * coef_grid * fft_2d(in(i, :, :)))
      end do
      !$OMP END PARALLEL DO
   end function partial_fft

   function curl(in, inv_sqrtg) result(out)
      use global, only : fth, fph, d_s => delta_s
      real(r8), intent(in) :: inv_sqrtg(:,:,:)
      type(vector_grid), intent(in):: in
      type(vector_grid) :: out, grad_1, grad_2, grad_3

      out%u1 = (partial_fft(in%u3, fth) - partial_fft(in%u2, fph)) * inv_sqrtg
      out%u2 = (partial_fft(in%u1, fph) - partial_rad(in%u3, d_s)) * inv_sqrtg
      out%u3 = (partial_rad(in%u2, d_s) - partial_fft(in%u1, fth)) * inv_sqrtg

   end function curl

   pure function finite_differences(y, dx) result(dy)
      complex(r8), intent(in) :: y(:)
      real(r8), intent(in) :: dx
      complex(r8) :: dy(size(y))
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
end module derivatives
