module types
   use constants

   implicit none

   type vector_grid
      real(8), allocatable, dimension(:,:,:) :: u1, u2, u3
   contains
      procedure :: init_grid => init_vector_grid
      procedure :: init_vec => init_vector_vectors
      procedure :: alloc => alloc_vector_grid
   end type vector_grid


   type metric_tensor
      real(8), allocatable, dimension(:,:,:) :: g11, g12, g13, g21, g22, g23, g31, g32, g33
   contains
      procedure :: init => init_metric_tensor
   end type metric_tensor

contains


   subroutine alloc_vector_grid(self, ls, lth, lph)
      class(vector_grid) :: self
      integer(8) :: ls, lth, lph

      allocate(self%u1(ls, lth, lph), &
         self%u2(ls, lth, lph), &
         self%u3(ls, lth, lph))

   end subroutine alloc_vector_grid

   subroutine init_vector_grid(self, grid)
      real(8), intent(in) :: grid(:,:,:,:)
      class(vector_grid) :: self

      self%u1 = grid(1,:,:,:)
      self%u2 = grid(2,:,:,:)
      self%u3 = grid(3,:,:,:)

   end subroutine init_vector_grid

   subroutine init_vector_vectors(self, v1, v2, v3)
      real(8), intent(in), dimension(:,:,:) :: v1, v2, v3
      class(vector_grid) :: self

      self%u1 = v1
      self%u2 = v2
      self%u3 = v3

   end subroutine init_vector_vectors

   elemental function cross_vector_grid(vec_a, vec_b) result(vec_c)
      type(vector_grid), intent(in) :: vec_a, vec_b
      type(vector_grid) :: vec_c

      vec_c%u1 = vec_a%u2 * vec_b%u3 - vec_a%u3 * vec_b%u2
      vec_c%u2 = vec_a%u3 * vec_b%u1 - vec_a%u1 * vec_b%u3
      vec_c%u3 = vec_a%u1 * vec_b%u2 - vec_a%u2 * vec_b%u1

   end function cross_vector_grid

   pure function sum_three_vector_grids(vec1, vec2, vec3) result(out)
      type(vector_grid), intent(in) :: vec1, vec2, vec3
      type(vector_grid) :: out

      out%u1 = vec1%u1 + vec2%u1 + vec3%u1
      out%u2 = vec1%u2 + vec2%u2 + vec3%u2
      out%u3 = vec1%u3 + vec2%u3 + vec3%u3
   end function sum_three_vector_grids

   pure function product_with_scalar_grid(vec_a, sc_grid) result(vec_b)
      type(vector_grid), intent(in) :: vec_a
      real(8), dimension(:,:,:), intent(in) :: sc_grid
      type(vector_grid):: vec_b

      vec_b%u1 = vec_a%u1 * sc_grid
      vec_b%u2 = vec_a%u2 * sc_grid
      vec_b%u3 = vec_a%u3 * sc_grid
   end function product_with_scalar_grid

   function product_with_scalar(vec_a, scalar) result(vec_b)
      type(vector_grid), intent(in) :: vec_a
      real(8), intent(in) :: scalar
      type(vector_grid):: vec_b

      vec_b%u1 = vec_a%u1 * scalar
      vec_b%u2 = vec_a%u2 * scalar
      vec_b%u3 = vec_a%u3 * scalar
   end function product_with_scalar

   function dot_vector_grid(vec_a, vec_b) result(grid_c)
      type(vector_grid), intent(in) :: vec_a, vec_b
      real(8) :: grid_c(size(vec_a%u1, 1), size(vec_a%u1, 2), size(vec_a%u1, 3))

      grid_c = vec_a%u1 * vec_b%u1 + vec_a%u2 * vec_b%u2 + vec_a%u3 * vec_b%u3

   end function dot_vector_grid

   subroutine init_metric_tensor(self, e_sub_s, e_sub_th, e_sub_ph)
      class(metric_tensor) :: self
      type(vector_grid) :: e_sub_s, e_sub_th, e_sub_ph

      self%g11 = dot_vector_grid(e_sub_s, e_sub_s)
      self%g12 = dot_vector_grid(e_sub_s, e_sub_th)
      self%g13 = dot_vector_grid(e_sub_s, e_sub_ph)

      self%g21 = dot_vector_grid(e_sub_th, e_sub_s)
      self%g22 = dot_vector_grid(e_sub_th, e_sub_th)
      self%g23 = dot_vector_grid(e_sub_th, e_sub_ph)

      self%g31 = dot_vector_grid(e_sub_ph, e_sub_s)
      self%g32 = dot_vector_grid(e_sub_ph, e_sub_th)
      self%g33 = dot_vector_grid(e_sub_ph, e_sub_ph)

   end subroutine init_metric_tensor

   function lower_indices(vec_super, g_ij) result(vec_sub)
      type(vector_grid), intent(in) :: vec_super
      type(metric_tensor), intent(in) :: g_ij
      type(vector_grid) :: vec_sub

      vec_sub%u1 = vec_super%u1 * g_ij%g11 + vec_super%u2 * g_ij%g12 + vec_super%u3 * g_ij%g13
      vec_sub%u2 = vec_super%u1 * g_ij%g21 + vec_super%u2 * g_ij%g22 + vec_super%u3 * g_ij%g23
      vec_sub%u3 = vec_super%u1 * g_ij%g31 + vec_super%u2 * g_ij%g32 + vec_super%u3 * g_ij%g33

   end function lower_indices

end module types
