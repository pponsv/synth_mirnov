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


   function partial_rad(in, ds) result(out)
      complex(r8), intent(in) :: in(:,:,:)
      real(r8) :: ds
      complex(r8) :: out(size(in, 1), size(in, 2), size(in, 3))
      integer :: j, k

      !$OMP PARALLEL DO PRIVATE(k, j)
      do k=1, size(in, 3)
         do j=1, size(in, 2)
            out(:, j, k) = finite_differences_1(in(:, j, k), ds)
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


   pure function finite_differences_1(y, dx) result(dy)
      complex(r8), intent(in) :: y(:)
      real(r8), intent(in) :: dx
      complex(r8) :: dy(size(y))
      integer :: i, len_y

      len_y = size(y)

      dy(1) = (y(2) - y(1)) / (dx)
      do i=2, len_y
         dy(i) = (y(i+1) - y(i-1)) / (2*dx)
      end do
      dy(len_y) = (y(len_y) - y(len_y - 1)) / (dx)

   end function finite_differences_1


   pure function finite_differences_2(y, dx) result(dy)
      complex(r8), intent(in) :: y(:)
      real(r8), intent(in) :: dx
      complex(r8) :: dy(size(y))
      integer :: i, len_y

      len_y = size(y)

      dy(1) = (-3*y(1) + 4*y(2) - 1*y(3)) / (2*dx)
      dy(2) = (-2*y(1) - 3*y(2) + 6*y(3) - 1*y(4)) / (6*dx)
      do i=3, len_y-2
         dy(i) = (1*y(i - 2) - 8*y(i - 1) + 8*y(i + 1) - 1*y(i + 2)) / (12*dx)
      end do
      dy(len_y - 1) = (1*y(len_y - 3) - 6*y(len_y - 2) + 3*y(len_y - 1) + 2*y(len_y)) / (6*dx)
      dy(len_y) = (1*y(len_y - 2) - 4*y(len_y - 1) + 3*y(len_y)) / (2*dx)

   end function finite_differences_2


   pure function finite_differences_3(y, dx) result(dy)
      complex(r8), intent(in) :: y(:)
      real(r8), intent(in) :: dx
      complex(r8) :: dy(size(y))
      integer :: i, len_y

      len_y = size(y)

      dy(1) = (-11*y(1) + 18*y(2) - 9*y(3) + 2*y(4)) / (6*dx)
      dy(2) = (-3*y(1) - 10*y(2) + 18*y(3) - 6*y(4) + y(5)) / (12*dx)
      dy(3) = (3*y(1) - 30*y(2) - 20*y(3) + 60*y(4) - 15*y(5) + 2*y(6)) / (60*dx)
      do i=4, len_y-3
         dy(i) = (-y(i - 3) + 9*y(i-2) - 45*y(i-1) + 45*y(i+1) - 9*y(i+2) + y(i+3)) / (60*dx)
      end do
      dy(len_y - 2) = (-2*y(len_y - 5) + 15*y(len_y - 4) - 60*y(len_y - 3) + 20*y(len_y-2) + 30*y(len_y - 1) - 3*y(len_y)) / (60*dx)
      dy(len_y - 1) = (-y(len_y - 4) + 6*y(len_y - 3) - 18*y(len_y - 2) + 10*y(len_y - 1) + 3*y(len_y)) / (12*dx)
      dy(len_y) = (-2*y(len_y - 3) + 9*y(len_y - 2) - 18*y(len_y - 1) + 11*y(len_y)) / (6*dx)

   end function finite_differences_3


   pure function finite_differences_4(y, dx) result(dy)
      complex(r8), intent(in) :: y(:)
      real(r8), intent(in) :: dx
      complex(r8) :: dy(size(y))
      integer :: i, len_y

      len_y = size(y)

      dy(1) = (-25*y(1) + 48*y(2) - 36*y(3) + 16*y(4) - 3*y(5)) / (12*dx)
      dy(2) = (-12*y(1) - 65*y(2) + 120*y(3) - 60*y(4) + 20*y(5) -3*y(6)) / (60*dx)
      dy(3) = (2*y(1) - 24*y(2) - 35*y(3) + 80*y(4) - 30*y(5) + 8*y(6) - y(7)) / (60*dx)
      dy(4) = (-4*y(1) + 42*y(2) - 252*y(3) - 105*y(4) + 420*y(5) - 126*y(6) + 28*y(7) - 3*y(8)) / (420*dx)
      do i=5, len_y-4
         dy(i) = (3*y(i-4) - 32*y(i - 3) + 168*y(i-2) - 672*y(i-1) + 672*y(i+1) - 168*y(i+2) + 32*y(i+3) - 3*y(i+4)) / (840*dx)
      end do
      dy(len_y - 3) =  (3*y(len_y - 7) - 28*y(len_y - 6) + 126*y(len_y - 5) - 420*y(len_y-4) + 105*y(len_y - 3) + &
         252*y(len_y - 2) - 42*y(len_y - 1) + 4*y(len_y)) / (420*dx)
      dy(len_y - 2) = (1*y(len_y - 6) - 8*y(len_y - 5) + 30*y(len_y - 4) - 80*y(len_y-3) + 35*y(len_y - 2) + 24*y(len_y - 1) - 2*y(len_y)) / (60*dx)
      dy(len_y - 1) = (3*y(len_y - 5) - 20*y(len_y - 4) + 60*y(len_y - 3) - 120*y(len_y - 2) + 65*y(len_y - 1) + 12*y(len_y)) / (60*dx)
      dy(len_y) = (3*y(len_y - 4) - 16*y(len_y - 3) + 36*y(len_y - 2) - 48*y(len_y - 1) + 25*y(len_y)) / (12*dx)

   end function finite_differences_4

end module derivatives
