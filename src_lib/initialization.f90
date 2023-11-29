module initialization
   use constants

   implicit none

contains

   subroutine initialize_boozer(s_b, th_b, ph_b, phi_b_g, iota_b, b_mod_b, sqrt_g_b, x, y, z)
      use fft_mod, only: fftfreqs, meshgrid, plan_ffts
      use global
      real(r8), intent(in) :: s_b(:), th_b(:), ph_b(:), phi_b_g, iota_b(:)
      real(r8), intent(in), dimension(:,:,:) :: b_mod_b, sqrt_g_b, x, y, z
      integer :: i

      !  Coordinate grids
      s = s_b
      th = th_b
      ph = ph_b

      len_s = size(s)
      len_th = size(th)
      len_ph = size(ph)

      !  We assume a uniform grid. For non-uniform grids, considerable
      !  modification is needed
      delta_s = s(2) - s(1)
      delta_th = th(2) - th(1)
      delta_ph = ph(2) - ph(1)
      int_factor = delta_s * delta_th * delta_ph

      !  FFTs for derivatives
      th_freqs = 2*PI*fftfreqs(len_th, delta_th)
      ph_freqs = 2*PI*fftfreqs(len_ph, delta_ph)
      allocate(fth(len_th, len_ph), fph(len_th, len_ph))
      call meshgrid(ph_freqs, th_freqs, fph, fth)
      call plan_ffts

      !  Profiles, magnetic field and jacobian
      iota = iota_b
      mod_b = b_mod_b
      inv_mod_b2 = 1. / b_mod_b**2
      sqrt_g = sqrt_g_b
      inv_sqrt_g = 1. / sqrt_g

      allocate(b_super(len_s, len_th, len_ph, 3))
      allocate(xyz_grid(len_s, len_th, len_ph, 3))

      !  phi_b_g is the derivative of the flux wrt the radial coordinate,
      !  equal to the flux at the LCFS if the radial coordinate is s
      b_super(:,:,:,1) = 0
      b_super(:,:,:,3) = -phi_b_g/(2*pi*sqrt_g)
      do i=1, len_s
         b_super(i, :, :, 2) = b_super(i, :, :, 3) * iota(i)
      end do

      !  Cartesian grid
      xyz_grid(:,:,:,1) = x
      xyz_grid(:,:,:,2) = y
      xyz_grid(:,:,:,3) = z
   end subroutine initialize_boozer


   subroutine initialize_basis(e_sub_s_b, e_sub_th_b, e_sub_ph_b)
      use global, only : e_sub_s, e_sub_th, e_sub_ph, g_sub_ij
      use helper, only : metric_tensor, checknans
      real(r8), intent(in), dimension(:,:,:,:) :: e_sub_s_b, e_sub_th_b, e_sub_ph_b

      e_sub_s = e_sub_s_b
      e_sub_th = e_sub_th_b
      e_sub_ph = e_sub_ph_b

      g_sub_ij = metric_tensor(e_sub_s, e_sub_th, e_sub_ph)

   end subroutine initialize_basis


   subroutine initialize_potentials(profiles, ms, ns, fs, time)
      use global
      integer(i8), intent(in) :: ms(:), ns(:)
      real(r8), intent(in) :: fs(:), time(:)
      complex(r8), intent(in) :: profiles(:,:)

      num_modes = size(ms)
      ms_pot = ms
      ns_pot = ns
      ws_pot = 2.*PI*fs
      t = time
      len_t = size(time)
      prof_nm = profiles

      call try_dealloc_pot

      allocate(dpot_dph(len_s, len_th, len_ph, num_modes))
      allocate(dpot_dth(len_s, len_th, len_ph, num_modes))
      allocate(gradpar_pot_super(len_s, len_th, len_ph, 3, num_modes))
      allocate(j_super(len_s, len_th, len_ph, 3, num_modes))
      allocate(db_all(len_s, len_th, len_ph, 3, num_modes))

   contains

      subroutine try_dealloc_pot
         use global, only : dpot_dph, dpot_dth, gradpar_pot_super, j_super, db_all

         if (allocated(dpot_dph)) deallocate(dpot_dph)
         if (allocated(dpot_dth)) deallocate(dpot_dth)
         if (allocated(gradpar_pot_super)) deallocate(gradpar_pot_super)
         if (allocated(j_super)) deallocate(j_super)
         if (allocated(db_all)) deallocate(db_all)

      end subroutine try_dealloc_pot

   end subroutine initialize_potentials


   subroutine initialize_coils(xyzs, ncoils)
      use global, only : coil_xyz, num_coils, len_s, len_th, len_ph, r_coil
      use helper, only : scalar_product_real, cross_product_real, dot_product_real
      use global, only : coil_xyz, xyz_grid, e_sub_s, e_sub_ph, e_sub_th,&
         us, uth, uph, r_coil, sqrt_g, r_3_sqrtg

      integer(i8), intent(in) :: ncoils
      real(r8), intent(in) :: xyzs(ncoils, 3)
      integer :: idx_coil

      write (*, '(/, A)', advance="no") "COIL INITIALIZATION: "

      coil_xyz = transpose(xyzs)
      num_coils = ncoils

      if (allocated(r_coil)) deallocate(r_coil)
      allocate(r_coil(len_s, len_th, len_ph, 3, num_coils))
      allocate(us(len_s, len_th, len_ph, 3, num_coils))
      allocate(uth(len_s, len_th, len_ph, 3, num_coils))
      allocate(uph(len_s, len_th, len_ph, 3, num_coils))

      do idx_coil=1, num_coils
         r_coil(:,:,:,1, idx_coil) = coil_xyz(1,idx_coil) - xyz_grid(:,:,:,1)
         r_coil(:,:,:,2, idx_coil) = coil_xyz(2,idx_coil) - xyz_grid(:,:,:,2)
         r_coil(:,:,:,3, idx_coil) = coil_xyz(3,idx_coil) - xyz_grid(:,:,:,3)

         r_3_sqrtg = sqrt_g * sqrt(dot_product_real(r_coil(:,:,:,:,idx_coil), &
            r_coil(:,:,:,:,idx_coil)))**(-3)

         us(:,:,:,:, idx_coil)  = scalar_product_real(cross_product_real(e_sub_s,  &
            r_coil(:,:,:,:,idx_coil)), r_3_sqrtg)
         uth(:,:,:,:, idx_coil) = scalar_product_real(cross_product_real(e_sub_th, &
            r_coil(:,:,:,:,idx_coil)), r_3_sqrtg)
         uph(:,:,:,:, idx_coil) = scalar_product_real(cross_product_real(e_sub_ph, &
            r_coil(:,:,:,:,idx_coil)), r_3_sqrtg)

      end do

      write (*, '(/, A)') "DONE"

   end subroutine initialize_coils

end module initialization
