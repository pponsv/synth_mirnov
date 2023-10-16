module derivatives

   use constants

   implicit none

contains

   function gradient(in) result(out)
      use global, only : len_s, len_ph, len_th, delta_s, fth, fph
      complex(r8), intent(in) :: in(len_s, len_th, len_ph)
      complex(r8):: out(len_s, len_th, len_ph, 3)

      out(:,:,:,1) = partial_rad(in, delta_s)
      out(:,:,:,2) = partial_fft(in, fth)
      out(:,:,:,3) = partial_fft(in, fph)

   end function gradient


   function partial_rad(in, ds) result(out)
      complex(r8), intent(in) :: in(:,:,:)
      real(r8) :: ds
      complex(r8) :: out(size(in, 1), size(in, 2), size(in, 3))
      integer :: j, k

      !$OMP PARALLEL DO PRIVATE(k, j)
      do k=1, size(in, 3)
         do j=1, size(in, 2)
            out(:, j, k) = finite_differences(in(:, j, k), ds)
         end do
      end do
      !$OMP END PARALLEL DO

   end function partial_rad


   function partial_fft(arr_in, coef_grid) result(out)
      use fft_mod, only : fft_2d, ifft_2d
      complex(r8), intent(in) :: arr_in(:,:,:)
      real(r8) :: coef_grid(:,:)
      complex(r8) :: out(size(arr_in, 1), size(arr_in, 2), size(arr_in, 3))
      integer :: i

      !$OMP PARALLEL DO PRIVATE(i)
      do i=1, size(arr_in, 1)
         out(i, :, :) = ifft_2d(imag * coef_grid * fft_2d(arr_in(i, :, :)))
      end do
      !$OMP END PARALLEL DO

   end function partial_fft


   function curl(vec_in, inv_sqrtg) result(out)
      use global, only : fth, fph, d_s => delta_s
      complex(r8), intent(in) :: vec_in(:,:,:,:)
      real(r8), intent(in) :: inv_sqrtg(:,:,:)
      complex(r8) :: out(size(vec_in, 1), size(vec_in, 2), size(vec_in, 3), 3)

      out(:,:,:,1) = (partial_fft(vec_in(:,:,:,3), fth) - &
         partial_fft(vec_in(:,:,:,2), fph)) * inv_sqrtg
      out(:,:,:,2) = (partial_fft(vec_in(:,:,:,1), fph) - &
         partial_rad(vec_in(:,:,:,3), d_s)) * inv_sqrtg
      out(:,:,:,3) = (partial_rad(vec_in(:,:,:,2), d_s) - &
         partial_fft(vec_in(:,:,:,1), fth)) * inv_sqrtg

   end function curl


   pure function finite_differences(y, dx) result(dy)
      complex(r8), intent(in) :: y(:)
      real(r8), intent(in) :: dx
      complex(r8) :: dy(size(y))
      integer :: i, len_y

      len_y = size(y)

      dy(1) = (-3*y(1) + 4*y(2) - 1*y(3)) / (2*dx)
      dy(2) = (-2*y(1) - 3*y(2) + 6*y(3) - 1*y(4)) / (6*dx)
      do i=3, size(y)-2
         dy(i) = (1*y(i - 2) - 8*y(i - 1) + 8*y(i + 1) - 1*y(i + 2)) / (12*dx)
      end do
      dy(size(y) - 1) = (1*y(len_y - 3) - 6*y(len_y - 2) + 3*y(len_y - 1) + 2*y(len_y)) / (6*dx)
      dy(size(y)) = (1*y(len_y - 2) - 4*y(len_y - 1) + 3*y(len_y)) / (2*dx)

   end function finite_differences

end module derivatives
