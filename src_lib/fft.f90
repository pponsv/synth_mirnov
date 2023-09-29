module fft_mod
   use constants
   use, intrinsic :: iso_c_binding
   implicit none
   include 'fftw3.f03'

contains

   subroutine plan_ffts
      use global, only : fft_plan_2d, ifft_plan_2d, len_th, len_ph
      complex(r8) :: out(len_th, len_ph), in(len_th, len_ph)

      print *, 'PLANS', fft_plan_2d, ifft_plan_2d

      fft_plan_2d = fftw_plan_dft_2d(int(len_ph, 4), int(len_th, 4), in, out, &
         FFTW_FORWARD, FFTW_ESTIMATE)
      ifft_plan_2d = fftw_plan_dft_2d(int(len_ph, 4), int(len_th, 4), out, in, &
         FFTW_BACKWARD, FFTW_ESTIMATE)

      print *, 'PLANS', fft_plan_2d, ifft_plan_2d

   end subroutine plan_ffts

   function fft_2d(in) result(out)
      use global, only : len_th, len_ph, fft_plan_2d
      complex(r8) :: out(len_th, len_ph), in(len_th, len_ph)

      print *, size(in, 1), size(in, 2)
      call fftw_execute_dft(fft_plan_2d, in, out)
      print *, out(1,1)
      ! stop
   end function fft_2d

   function ifft_2d(in) result(out)
      use global, only : len_th, len_ph, ifft_plan_2d
      complex(r8) :: in(len_th, len_ph), out(len_th, len_ph)

      call fftw_execute_dft(ifft_plan_2d, in, out)
      out = out / (len_th * len_ph)

   end function ifft_2d

   function fftfreqs(n, d) result(freqs)
      integer(i8), intent(in) :: n
      real(r8), intent(in) :: d
      real(r8) :: freqs(n), fac
      integer(i8) :: i, n_half

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
      real(r8), intent(in) :: a(:), b(:)
      real(r8), intent(out) :: aa(size(b), size(a)), bb(size(b), size(a))
      integer :: i, j

      do i=1, size(b)
         aa(i, :) = a
      end do
      do i=1, size(a)
         bb(:, i) = b
      end do
   end subroutine meshgrid

end module fft_mod
