module helper

   use constants

   implicit none

contains

   subroutine test
      print *, "Test successful"
   end subroutine test

   function real_linspace(xi, xf, np) result(out)
      !       Makes a linearly spaced vector between xi, xf
      !       with np elements

      !   Input
      real(8), intent(in)    :: xi, xf
      integer(8), intent(in) :: np
      !   Output
      real(8)    :: out(np)
      !   Internal
      integer(8) :: i

      do i = 1, np
         out(i) = xi + (xf - xi)/(np - 1)*(i - 1)
      end do
   end function real_linspace

   subroutine genvar_new(fnm, th, ph, ms, ns, n_s, n_m, n_th, n_ph, typ, out)
      integer(8), intent(in) :: n_m, n_th, n_ph, n_s
      real(8), intent(in) :: fnm(n_m, n_s), th(n_th), ph(n_ph), ms(n_m), ns(n_m)
      character, intent(in) :: typ
      real(8), intent(out) :: out(n_s, n_th, n_ph)

      complex(8) :: tmp_out(n_s, n_th, n_ph)
      integer(8) :: i, j, k

      !$OMP PARALLEL DO PRIVATE(i, j, k)
      do k = 1, n_ph
         do j = 1, n_th
            do i = 1, n_s
               ! tmp_out(i, j, k) = sum(fnm(:, i)*exp_mth(:,j)*exp_nph(:,k))
            end do
         end do
      end do
      !$OMP END PARALLEL DO

      if (typ == 'c') then
         out = tmp_out % re
      else if (typ == 's') then
         out = tmp_out % im
      end if

   end subroutine genvar_new



   function dfnm_ds(fnm, s, n_m, n_s) result(dfnm)
      real(8), intent(in) :: s(n_s), fnm(n_m, n_s)
      integer(8), intent(in) :: n_m, n_s
      real(8) :: dfnm(n_m, n_s)
      integer(8) :: i
      !  Finite differences coefficients from:
      !  https://web.media.mit.edu/~crtaylor/calculator.html
      !$OMP PARALLEL DO PRIVATE(i)
      do i = 1, n_s
         if (i == 1) then
            ! dfnm(:, i) = (fnm(:, i + 1) - fnm(:, i))/((s(i + 1) - s(i))) ! First order forward
            dfnm(:, i) = (4*fnm(:, i + 1) - 3*fnm(:, i) - fnm(:, i + 2))/(2*(s(i + 1) - s(i))) ! Second order forward
         else if (i == n_s) then
            ! dfnm(:, i) = (fnm(:, i) - fnm(:, i - 1))/((s(i) - s(i - 1))) ! First order backwards
            dfnm(:, i) = (fnm(:,i-2) + 3*fnm(:, i) - 4*fnm(:, i - 1))/(2*(s(i) - s(i - 1))) ! Second order backwards
         else
            dfnm(:, i) = (fnm(:, i + 1) - fnm(:, i - 1))/((s(i + 1) - s(i - 1)))
         end if
      end do
      !$OMP END PARALLEL DO
   end function dfnm_ds



end module helper
