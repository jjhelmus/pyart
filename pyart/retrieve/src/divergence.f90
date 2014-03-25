!  Module: divergence.f90


subroutine horiz_wind(u, v, dx, dy, finite_scheme, fill_value, proc, &
                      nx, ny, nz, div, du, dv)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: dx, dy, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v
   real(kind=8), intent(out), dimension(nz,ny,nx) :: div, du, dv


!  Define local variables ====================================================

   real(kind=8), parameter :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                              b1=3.d0/4.d0, b2=3.d0/20.d0, &
                              b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                              c1=6.d0, c2=-15.d0/2.d0, &
                              c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                              c5=6.d0/5.d0, c6=-1.d0/6.d0

   integer(kind=4)         :: i, j, k

!  ===========================================================================


!  F2PY directives ============================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: u, v, dx, dy, fill_value
   !f2py real(kind=8), intent(out)             :: div, du, dv

!  ===========================================================================


!  The horizontal wind divergence is given by the sum of the 2 terms
!  du/dx and dv/dy, therefore we need to compute these 2 terms
!
!  We compute these partial derivatives in the so-called grid space (3-D)
!  rather than in the vector space because the problem is most natural in the
!  grid space
!
!  The first block is for the basic finite difference schemes

   !$omp parallel num_threads(proc)

   if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           The following is very important as it describes how we will
!           be approximating du/dx and dv/dy in order to compute the
!           horizontal wind divergence
!
!           For the interior points of the grid,
!
!           i = [2, nx-1] or
!           j = [2, ny-1]
!
!           use a centered difference scheme with p = 2. When at the
!           boundaries of the grid,
!
!           i = 1, nx or
!           j = 1, ny
!
!           use either a forward or backward difference scheme, both with
!           p = 1
!
!           First compute du/dx

            if (i > 1 .and. i < nx) then
               du(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

            elseif (i == 1) then
               du(k,j,i) = (u(k,j,i+1) - u(k,j,i)) / dx

            else
               du(k,j,i) = (u(k,j,i) - u(k,j,i-1)) / dx
            endif

!           Now compute dv/dy

            if (j > 1 .and. j < ny) then
               dv(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dv(k,j,i) = (v(k,j+1,i) - v(k,j,i)) / dy

            else
               dv(k,j,i) = (v(k,j,i) - v(k,j-1,i)) / dy
            endif

            enddo
         enddo
      enddo
      !$omp end do


!  The second block is for the high-order finite difference schemes

   elseif (finite_scheme == 'high-order') then
   
      !$omp do
      do i =  1, nx
         do j = 1, ny
            do k = 1, nz

!           The following is very important as it describes how we will
!           be approximating du/dx and dv/dy in order to compute the
!           horizontal wind divergence
!
!           For the interior points of the grid,
!
!           i = [4, nx-3] or
!           j = [4, ny-3]
!
!           use a centered difference scheme with p = 6. Closer to the
!           boundaries of the grid,
!
!           i = 2, 3, nx-2, nx-1 or
!           j = 2, 3, ny-2, ny-1
!
!           still use centered difference schemes, but of lower accuracy,
!           i.e. p = 2 or p = 4. When at the boundaries of the grid,
!
!           i = 1, nx or
!           j = 1, ny
!
!           use either a forward or backward difference scheme, both with
!           p = 6
!
!           First compute du/dx

            if (i > 3 .and. i < nx - 2) then
               du(k,j,i) = (b3 * u(k,j,i-3) + b2 * u(k,j,i-2) - &
                            b1 * u(k,j,i-1) + b1 * u(k,j,i+1) - &
                            b2 * u(k,j,i+2) - b3 * u(k,j,i+3)) / dx

            elseif (i > 2 .and. i < nx - 1) then
               du(k,j,i) = (a2 * u(k,j,i-2) - a1 * u(k,j,i-1) + &
                            a1 * u(k,j,i+1) - a2 * u(k,j,i+2)) / dx
                            
            elseif (i > 1 .and. i < nx) then
               du(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)
               
            elseif (i == 1) then
               du(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i+1) + &
                            c2 * u(k,j,i+2) + c3 * u(k,j,i+3) + &
                            c4 * u(k,j,i+4) + c5 * u(k,j,i+5) + &
                            c6 * u(k,j,i+6)) / dx
                            
            else
               du(k,j,i) = (-c0 * u(k,j,i) - c1 * u(k,j,i-1) - &
                             c2 * u(k,j,i-2) - c3 * u(k,j,i-3) - &
                             c4 * u(k,j,i-4) - c5 * u(k,j,i-5) - &
                             c6 * u(k,j,i-6)) / dx
            endif
            
!           Now compute dv/dy

            if (j > 3 .and. j < ny - 2) then
               dv(k,j,i) = (b3 * v(k,j-3,i) + b2 * v(k,j-2,i) - &
                            b1 * v(k,j-1,i) + b1 * v(k,j+1,i) - &
                            b2 * v(k,j+2,i) - b3 * v(k,j+3,i)) / dy
            
            elseif (j > 2 .and. j < ny - 1) then
               dv(k,j,i) = (a2 * v(k,j-2,i) - a1 * v(k,j-1,i) + &
                            a1 * v(k,j+1,i) - a2 * v(k,j+2,i)) / dy
                            
            elseif (j > 1 .and. j < ny) then
               dv(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)
               
            elseif (j == 1) then
               dv(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j+1,i) + &
                            c2 * v(k,j+2,i) + c3 * v(k,j+3,i) + &
                            c4 * v(k,j+4,i) + c5 * v(k,j+5,i) + &
                            c6 * v(k,j+6,i)) / dy

            else
               dv(k,j,i) = (-c0 * v(k,j,i) - c1 * v(k,j-1,i) - &
                             c2 * v(k,j-2,i) - c3 * v(k,j-3,i) - &
                             c4 * v(k,j-4,i) - c5 * v(k,j-5,i) - &
                             c6 * v(k,j-6,i)) / dy
            endif
            
            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif
   
!  Add du/dx and dv/dy to give the horizontal wind divergence

   div = du + dv

   !$omp end parallel

   return

end subroutine horiz_wind


subroutine full_wind(u, v, w, dx, dy, dz, finite_scheme, fill_value, proc, &
                     nx, ny, nz, div, du, dv, dw)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w
   real(kind=8), intent(out), dimension(nz,ny,nx) :: div, du, dv, dw

!  Local variables ===========================================================

   real(kind=8), parameter :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                              b1=3.d0/4.d0, b2=3.d0/20.d0, &
                              b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                              c1=6.d0, c2=-15.d0/2.d0, &
                              c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                              c5=6.d0/5.d0, c6=-1.d0/6.d0

   integer(kind=4)         :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w
   !f2py real(kind=8), intent(out)             :: div, du, dv, dw

!  ===========================================================================


!  The 3-D divergence is given by the sum of the 3 terms du/dx, dv/dy
!  and dw/dz. Therefore we need to compute these 3 terms
!
!  We compute these partial derivatives in the so-called grid space (3-D)
!  rather than in the vector space because the problem is most natural in the
!  grid space
!
!  The first block is for the basic finite difference schemes

   !$omp parallel num_threads(proc)

   if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           The following is very important as it describes how we will
!           be approximating du/dx, dv/dy, and dw/dz in order to compute the
!           3-D wind divergence
!
!           For the interior points of the grid,
!
!           i = [2, nx-1] or
!           j = [2, ny-1] or
!           k = [2, nz-1]
!
!           use a centered difference scheme with p = 2. When at the
!           boundaries of the grid,
!
!           i = 1, nx or
!           j = 1, ny or
!           k = 1, nz
!
!           use either a forward or backward difference scheme, both with
!           p = 1
!
!           First compute du/dx

            if (i > 1 .and. i < nx) then
               du(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

             elseif (i == 1) then
               du(k,j,i) = (u(k,j,i+1) - u(k,j,i)) / dx

            else
               du(k,j,i) = (u(k,j,i) - u(k,j,i-1)) / dx

            endif

!           Now compute dv/dy

            if (j > 1 .and. j < ny) then
               dv(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dv(k,j,i) = (v(k,j+1,i) - v(k,j,i)) / dy

            else
               dv(k,j,i) = (v(k,j,i) - v(k,j-1,i)) / dy
            endif

!           Now compute dw/dz

            if (k > 1 .and. k < nz) then
               dw(k,j,i) = (w(k+1,j,i) - w(k-1,j,i)) / (2.d0 * dz)

            elseif (k == 1) then
               dw(k,j,i) = (w(k+1,j,i) - w(k,j,i)) / dz

            else
               dw(k,j,i) = (w(k,j,i) - w(k-1,j,i)) / dz
            endif

            enddo
         enddo
      enddo
      !$omp end do
   

!  The second block is for the high-order finite difference schemes

   elseif (finite_scheme == 'high-order') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nx

!           The following is very important as it describes how we will
!           be approximating du/dx, dv/dy, and dw/dz in order to compute the
!           3-D wind divergence
!
!           For the interior points of the grid,
!
!           i = [4, nx-3] or
!           j = [4, ny-3] or
!           k = [4, nz-3]
!
!           use a centered difference scheme with p = 6. When closer to the
!           boundaries of the grid,
!
!           i = 2, 3, nx-2, nx-1 or
!           j = 2, 3, ny-2, ny-1 or
!           k = 2, 3, nz-2, nz-1
!
!           still use centered difference schemes, but of lower accuracy,
!           i.e. p = 2 or p = 4. When at the boundaries of the grid,
!
!           i = 1, nx or
!           j = 1, ny or
!           k = 1, nz
!
!           use either a forward or backward difference scheme, both with
!           p = 6
!
!           First compute du/dx

            if (i > 3 .and. i < nx - 2) then
               du(k,j,i) = (b3 * u(k,j,i-3) + b2 * u(k,j,i-2) - &
                            b1 * u(k,j,i-1) + b1 * u(k,j,i+1) - &
                            b2 * u(k,j,i+2) - b3 * u(k,j,i+3)) / dx

            elseif (i > 2 .and. i < nx - 1) then
               du(k,j,i) = (a2 * u(k,j,i-2) - a1 * u(k,j,i-1) + &
                            a1 * u(k,j,i+1) - a2 * u(k,j,i+2)) / dx

            elseif (i > 1 .and. i < nx) then
               du(k,j,i) = (u(k,j,i+1) - u(k,j,i-1)) / (2.d0 * dx)

            elseif (i == 1) then
               du(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i+1) + &
                            c2 * u(k,j,i+2) + c3 * u(k,j,i+3) + &
                            c4 * u(k,j,i+4) + c5 * u(k,j,i+5) + &
                            c6 * u(k,j,i+6)) / dx

            else
               du(k,j,i) = (-c0 * u(k,j,i) - c1 * u(k,j,i-1) - &
                             c2 * u(k,j,i-2) - c3 * u(k,j,i-3) - &
                             c4 * u(k,j,i-4) - c5 * u(k,j,i-5) - &
                             c6 * u(k,j,i-6)) / dx
            endif

!           Now compute dv/dy

            if (j > 3 .and. j < ny - 2) then
               dv(k,j,i) = (b3 * v(k,j-3,i) + b2 * v(k,j-2,i) - &
                            b1 * v(k,j-1,i) + b1 * v(k,j+1,i) - &
                            b2 * v(k,j+2,i) - b3 * v(k,j+3,i)) / dy

            elseif (j > 2 .and. j < ny - 1) then
               dv(k,j,i) = (a2 * v(k,j-2,i) - a1 * v(k,j-1,i) + &
                            a1 * v(k,j+1,i) - a2 * v(k,j+2,i)) / dy

            elseif (j > 1 .and. j < ny) then
               dv(k,j,i) = (v(k,j+1,i) - v(k,j-1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dv(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j+1,i) + &
                            c2 * v(k,j+2,i) + c3 * v(k,j+3,i) + &
                            c4 * v(k,j+4,i) + c5 * v(k,j+5,i) + &
                            c6 * v(k,j+6,i)) / dy

            else
               dv(k,j,i) = (-c0 * v(k,j,i) - c1 * v(k,j-1,i) - &
                             c2 * v(k,j-2,i) - c3 * v(k,j-3,i) - &
                             c4 * v(k,j-4,i) - c5 * v(k,j-5,i) - &
                             c6 * v(k,j-6,i)) / dy
            endif

!           Now compute dw/dz

            if (k > 3 .and. k < nz - 2) then
               dw(k,j,i) = (b3 * w(k-3,j,i) + b2 * w(k-2,j,i) - &
                            b1 * w(k-1,j,i) + b1 * w(k+1,j,i) - &
                            b2 * w(k+2,j,i) - b3 * w(k+3,j,i)) / dz

            elseif (k > 2 .and. k < nz - 1) then
               dw(k,j,i) = (a2 * w(k-2,j,i) - a1 * w(k-1,j,i) + &
                            a1 * w(k+1,j,i) - a2 * w(k+2,j,i)) / dz

            elseif (k > 1 .and. k < nz) then
               dw(k,j,i) = (w(k+1,j,i) - w(k-1,j,i)) / (2.d0 * dz)

            elseif (k == 1) then
               dw(k,j,i) = (c0 * w(k,j,i) + c1 * w(k+1,j,i) + &
                            c2 * w(k+2,j,i) + c3 * w(k+3,j,i) + &
                            c4 * w(k+4,j,i) + c5 * w(k+5,j,i) + &
                            c6 * w(k+6,j,i)) / dz

            else
               dw(k,j,i) = (-c0 * w(k,j,i) - c1 * w(k-1,j,i) - &
                             c2 * w(k-2,j,i) - c3 * w(k-3,j,i) - &
                             c4 * w(k-4,j,i) - c5 * w(k-5,j,i) - &
                             c6 * w(k-6,j,i)) / dz
            endif

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop
   
   endif
   
!  Add du/dx, dv/dy, and dw/dz to give the 3-D wind divergence

   div = du + dv + dw

   !$omp end parallel

   return

end subroutine full_wind


subroutine sub_beam(div, base, column, z, fill_value, proc, nx, ny, nz)

   implicit none

   integer(kind=4), intent(in)                      :: nx, ny, nz, proc
   real(kind=8), intent(in)                         :: fill_value
   real(kind=8), intent(in), dimension(nz)          :: z
   real(kind=8), intent(in), dimension(ny,nx)       :: base
   integer(kind=4), intent(in), dimension(ny,nx)    :: column
   real(kind=8), intent(inout), dimension(nz,ny,nx) :: div


!  Local variables ===========================================================

   integer(kind=4)           :: i, j, k

!  ============================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: column, proc
   !f2py real(kind=8), intent(in)              :: base, z, fill_value
   !f2py real(kind=8), intent(in,out)          :: div

!  ============================================================================

   
!  Implement sub-beam divergence criteria
!
!  Wherever the column type is well-defined or top-defined, he horizontal
!  wind divergence of each grid point near the surface with missing data is
!  set equal to the horizontal wind divergence of the first valid grid point
!  in the same column
!
!  The types of columns are as follows:
!
!  0 = Undefined
!  1 = Well-defined
!  2 = Top-defined
!  3 = Anvil-like
!  4 = Transition-like
!  5 = Discontinuous

   !$omp parallel num_threads(proc)

   !$omp do private(k)
   do i = 1, nx
      do j = 1, ny

!        Column types 1 and 2 are reserved for well-defined and top-defined,
!        respectively. If the current column is flagged with either of these
!        values, then we look to see where the echo base is, and set the
!        horizontal wind divergence below this height equal to that at the
!        echo base

         if (column(j,i) == 1 .or. column(j,i) == 2) then

         if (base(j,i) > z(1)) then
            k = minloc(abs(z - base(j,i)), dim=1)
            div(1:k,j,i) = div(k,j,i)
         endif

         endif

      enddo
   enddo
   !$omp end do
   
   !$omp end parallel

   return

end subroutine sub_beam
