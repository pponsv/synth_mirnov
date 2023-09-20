module fft_mod
   use, intrinsic :: iso_c_binding
   implicit none
   include 'fftw3.f03'

contains

   subroutine plan_ffts()
      use global, only : fft_plan_2d, ifft_plan_2d, len_th, len_ph

      real(8) :: in(len_th, len_ph)
      complex(8) :: out(len_th, len_ph)

      fft_plan_2d = fftw_plan_dft_r2c_2d(int(len_ph, 4), int(len_th, 4), in, out, &
         FFTW_ESTIMATE)
      ifft_plan_2d = fftw_plan_dft_c2r_2d(int(len_ph, 4), int(len_th, 4), out, in, &
         FFTW_ESTIMATE)

   end subroutine plan_ffts

   function fft_2d(in) result(out)
      use global, only : len_th, len_ph, fft_plan_2d
      real(8) :: in(len_th, len_ph)
      complex(8) :: out(len_th, len_ph)

      call fftw_execute_dft_r2c(fft_plan_2d, in, out)

   end function fft_2d

   function ifft_2d(in) result(out)
      use global, only : len_th, len_ph, ifft_plan_2d
      complex(8) :: in(len_th, len_ph)
      real(8) :: out(len_th, len_ph)

      call fftw_execute_dft_c2r(ifft_plan_2d, in, out)
      out = out / (len_th * len_ph)

   end function ifft_2d

   function fftfreqs(n, d) result(freqs)
      integer(8), intent(in) :: n
      real(8), intent(in) :: d
      real(8) :: freqs(n), fac
      integer(8) :: n_half
      integer(8) :: i

      fac = 1.0 / (n * d)

      n_half = (n-1)/2 +1
      do i=1, n_half
         freqs(i) = i-1
      end do
      do i=1, n-n_half
         freqs(i+n_half) = -(n/2) + i -1
      end do

      freqs = freqs * fac

   end function fftfreqs

   subroutine meshgrid(a, b, aa, bb)
      real(8), intent(in) :: a(:), b(:)
      real(8), intent(out) :: aa(size(a), size(b)), bb(size(a), size(b))
      integer :: i, j

      do i=1, size(b)
         aa(:, i) = a
      end do
      do i=1, size(a)
         bb(i, :) = b
      end do
   end subroutine meshgrid

end module fft_mod
