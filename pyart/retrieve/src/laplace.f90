!  Module: laplace.f90


subroutine orig_wind(u, v, dx, dy, dz, finite_scheme, fill_value, proc, &
                     ny, nx, nz, dux, duy, duz, duxy, duxz, duyz, &
                     dvx, dvy, dvz, dvxy, dvxz, dvyz)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dux, duy, duz, duxy, &
                                                     duxz, duyz, dvx, dvy, &
                                                     dvz, dvxy, dvxz, dvyz

!  Define local variables ====================================================

   real(kind=8)    :: dxx, dyy, dzz, dxy, dxz, dyz

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
   !f2py real(kind=8), intent(in)              :: u, v
   !f2py real(kind=8), intent(out)             :: dux, duy, duw, duxy, duxz
   !f2py real(kind=8), intent(out)             :: duyz, dvx, dvy, dvw, dvxy
   !f2py real(kind=8), intent(out)             :: dvxz, dvyz

!  ===========================================================================


!  Here we compute the vector Laplacian of the horizontal wind components
!  u and v, which means we have a total of 6 terms to compute at each grid
!  point. These 6 terms are d2u/dx2, d2u/dy2, d2u/dz2, d2v/dx2, d2v/dy2,
!  and d2v/dz2
!
!  We also compute the mixed derivatives of the horizontal wind components
!  (u,v), which means we have another 6 terms to compute at each grid point.
!  These 6 terms are d2u/dxdy, d2u/dxdz, d2u/dydz, d2v/dxdy, d2v/dxdz,
!  and d2v/dydz
!
!  We compute these partial derivatives in the so-called grid space (3-D)
!  rather than the vector space because it is most natural in the grid
!  space


!  First compute the parameters that will be used in the finite
!  differences

   dxx = dx**2
   dyy = dy**2
   dzz = dz**2

   dxy = dx * dy
   dxz = dx * dz
   dyz = dy * dz


!  The first block is for basic finite difference schemes

   !$omp parallel num_threads(proc)

   if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           First we will compute all the terms related to the vector
!           Laplacian
!
!           The following is very important as it describes how we will
!           be approximating the 6 vector Laplacian terms using finite
!           differences
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
!           First compute all the terms involving x derivatives, so here we
!           calculate d2u/dx2 and d2v/dx2

            if (i > 1 .and. i < nx) then
               dux(k,j,i) = (u(k,j,i+1) - 2.d0 * u(k,j,i) + u(k,j,i-1)) / dxx
               dvx(k,j,i) = (v(k,j,i+1) - 2.d0 * v(k,j,i) + v(k,j,i-1)) / dxx

            elseif (i == 1) then
               dux(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j,i+1) + u(k,j,i+2)) / dxx
               dvx(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j,i+1) + v(k,j,i+2)) / dxx

            else
               dux(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j,i-1) + u(k,j,i-2)) / dxx
               dvx(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j,i-1) + v(k,j,i-2)) / dxx
            endif

!           Now compute all the terms involving y derivatives, so here we
!           calculate d2u/dy2 and d2v/dy2

            if (j > 1 .and. j < ny) then
               duy(k,j,i) = (u(k,j+1,i) - 2.d0 * u(k,j,i) + u(k,j-1,i)) / dyy
               dvy(k,j,i) = (v(k,j+1,i) - 2.d0 * v(k,j,i) + v(k,j-1,i)) / dyy

            elseif (j == 1) then
               duy(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j+1,i) + u(k,j+2,i)) / dyy
               dvy(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j+1,i) + v(k,j+2,i)) / dyy

            else
               duy(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j-1,i) + u(k,j-2,i)) / dyy
               dvy(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j-1,i) + v(k,j-2,i)) / dyy
            endif

!           Now compute all the terms involving z derivatives, so here we
!           calculate d2u/dz2 and d2v/dz2

            if (k > 1 .and. k < nz) then
               duz(k,j,i) = (u(k+1,j,i) - 2.d0 * u(k,j,i) + u(k-1,j,i)) / dzz
               dvz(k,j,i) = (v(k+1,j,i) - 2.d0 * v(k,j,i) + v(k-1,j,i)) / dzz

            elseif (k == 1) then
               duz(k,j,i) = (u(k,j,i) - 2.d0 * u(k+1,j,i) + u(k+2,j,i)) / dzz
               dvz(k,j,i) = (v(k,j,i) - 2.d0 * v(k+1,j,i) + v(k+2,j,i)) / dzz

            else
               duz(k,j,i) = (u(k,j,i) - 2.d0 * u(k-1,j,i) + u(k-2,j,i)) / dzz
               dvz(k,j,i) = (v(k,j,i) - 2.d0 * v(k-1,j,i) + v(k-2,j,i)) / dzz
            endif

!           Now we want to compute all the terms related to the mixed
!           derivatives. Mixed derivatives require more care with respect
!           to grid boundaries since finite difference stencils have to
!           be used in multiple dimensions for every grid point
!
!           First compute all the terms involving x and y derivatives, so
!           here we calculate d2u/dxdy and d2u/dxdy
!
!           For the interior grid points for both x and y,
!
!           i = [2, nx-1] and
!           j = [2, ny-1]
!
!           use a centered difference scheme for both x and y derivatives,
!           with p = 2 for both. When the grid point is an interior point
!           for x but a boundary point for y,
!
!           i = [2, nx-1] and
!           j = 1, ny
!
!           then use a centered difference scheme for x with p = 2 and a
!           forward or backward scheme for y with p = 1. Use a similar
!           approach when grid point is interior for y but at the boundary
!           for x. Finally, when at the boundaries of the grid for both x
!           and y,
!
!           i = 1, nx and
!           j = 1, ny
!
!           use either a forward or backward difference scheme, both with
!           p = 1

            if (j > 1 .and. j < ny - 1 .and. i > 1 .and. i < nx - 1) then
               duxy(k,j,i) = (u(k,j+1,i+1) - u(k,j-1,i+1) - u(k,j+1,i-1) + &
                              u(k,j-1,i-1)) / (4.d0 * dxy)
               dvxy(k,j,i) = (v(k,j+1,i+1) - v(k,j-1,i+1) - v(k,j+1,i-1) + &
                              v(k,j-1,i-1)) / (4.d0 * dxy)

            elseif (j > 1 .and. j < ny - 1 .and. i == 1) then
               duxy(k,j,i) = (u(k,j+1,i+1) - u(k,j-1,i+1) - u(k,j+1,i) + &
                              u(k,j-1,i)) / (2.d0 * dxy)
               dvxy(k,j,i) = (v(k,j+1,i+1) - v(k,j-1,i+1) - v(k,j+1,i) + &
                              v(k,j-1,i)) / (2.d0 * dxy)

            elseif (j > 1 .and. j < ny - 1 .and. i == nx) then
               duxy(k,j,i) = (u(k,j+1,i) - u(k,j-1,i) - u(k,j+1,i-1) + &
                              u(k,j-1,i-1)) / (2.d0 * dxy)
               dvxy(k,j,i) = (v(k,j+1,i) - v(k,j-1,i) - v(k,j+1,i-1) + &
                              v(k,j-1,i-1)) / (2.d0 * dxy)

            elseif (j == 1 .and. i > 1 .and. i < nx - 1) then
               duxy(k,j,i) = (u(k,j+1,i+1) - u(k,j,i+1) - u(k,j+1,i-1) + &
                              u(k,j,i-1)) / (2.d0 * dxy)
               dvxy(k,j,i) = (v(k,j+1,i+1) - v(k,j,i+1) - v(k,j+1,i-1) + &
                              v(k,j,i-1)) / (2.d0 * dxy)

            elseif (j == ny .and. i > 1 .and. i < nx - 1) then
               duxy(k,j,i) = (u(k,j,i+1) - u(k,j-1,i+1) - u(k,j,i-1) + &
                              u(k,j-1,i-1)) / (2.d0 * dxy)
               dvxy(k,j,i) = (v(k,j,i+1) - v(k,j-1,i+1) - v(k,j,i-1) + &
                              v(k,j-1,i-1)) / (2.d0 * dxy)

            elseif (j == 1 .and. i == 1) then
               duxy(k,j,i) = (u(k,j,i) - u(k,j,i+1) - u(k,j+1,i) + &
                              u(k,j+1,i+1)) / dxy
               dvxy(k,j,i) = (v(k,j,i) - v(k,j,i+1) - v(k,j+1,i) + &
                              v(k,j+1,i+1)) / dxy

            else
               duxy(k,j,i) = (u(k,j,i) - u(k,j-1,i) - u(k,j,i-1) + &
                              u(k,j-1,i-1)) / dxy
               dvxy(k,j,i) = (v(k,j,i) - v(k,j-1,i) - v(k,j,i-1) + &
                              v(k,j-1,i-1)) / dxy
            endif

            enddo
         enddo
      enddo
      !$omp end do


!  The second block is for high-order finite difference schemes

   elseif (finite_scheme == 'high-order') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nx

!           TODO: Add capabilities for using high-order finite difference
!           schemes to approximate the partial derivatives

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

   return

end subroutine orig_wind


subroutine full_wind(u, v, w, dx, dy, dz, finite_scheme, fill_value, proc, &
                     nx, ny, nz, dux, duy, duz, dvx, dvy, dvz, dwx, &
                     dwy, dwz)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dux, duy, duz, dvx, &
                                                     dvy, dvz, dwx, dwy, &
                                                     dwz


!  Define local variables ====================================================

   real(kind=8), parameter :: a0=-5.d0/4.d0, a1=2.d0/3.d0, &
                              a2=-1.d0/24.d0, b0=-49.d0/36.d0, &
                              b1=3.d0/4.d0, b2=-3.d0/40.d0, &
                              b3=1.d0/180.d0, c0=469.d0/180.d0, &
                              c1=-223.d0/20.d0, c2=879.d0/40.d0, &
                              c3=-949.d0/36.d0, c4=41.d0/2.d0, &
                              c5=-201.d0/20.d0, c6=1019.d0/360.d0, &
                              c7=-7.d0/20.d0

   real(kind=8)            :: dxx, dyy, dzz

   integer(kind=4)         :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w
   !f2py real(kind=8), intent(out)             :: dux, duy, duz, dvx, dvy
   !f2py real(kind=8), intent(out)             :: dvz, dwx, dwy, dwz

!  ===========================================================================


!  Here we compute the vector Laplacian of the 3 wind components u, v, and w,
!  which means we have a total of 9 terms do compute at each grid point.
!  These 9 terms are d2u/dx2, d2u/dy2, d2u/dz2, d2v/dx2, d2v/dy2, d2v/dz2,
!  d2w/dx2, d2w/dy2, and d2w/dz2
!
!  We compute these partial derivatives in the so-called grid space (3-D)
!  rather than the vector space because it is most natural in the grid
!  space

!  First compute the parameters that will be used in the finite
!  differences

   dxx = dx**2
   dyy = dy**2
   dzz = dz**2

   !$omp parallel num_threads(proc)

!  The first block is for basic finite difference schemes

   if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           The following is very important as it describes how we will
!           be approximating the 9 vector Laplacian terms using finite
!           differences
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
!           First compute all the terms involving x derivatives, so here we
!           calculate d2u/dx2 and d2v/dx2, and d2w/dx2

            if (i > 1 .and. i < nx) then
               dux(k,j,i) = (u(k,j,i+1) - 2.d0 * u(k,j,i) + u(k,j,i-1)) / dxx
               dvx(k,j,i) = (v(k,j,i+1) - 2.d0 * v(k,j,i) + v(k,j,i-1)) / dxx
               dwx(k,j,i) = (w(k,j,i+1) - 2.d0 * w(k,j,i) + w(k,j,i-1)) / dxx

            elseif (i == 1) then
               dux(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j,i+1) + u(k,j,i+2)) / dxx
               dvx(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j,i+1) + v(k,j,i+2)) / dxx
               dwx(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j,i+1) + w(k,j,i+2)) / dxx

            else
               dux(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j,i-1) + u(k,j,i-2)) / dxx
               dvx(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j,i-1) + v(k,j,i-2)) / dxx
               dwx(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j,i-1) + w(k,j,i-2)) / dxx
            endif

!           Now compute all the terms involving y derivatives, so here we
!           calculate d2u/dy2, d2v/dy2, and d2w/dy2

            if (j > 1 .and. j < ny) then
               duy(k,j,i) = (u(k,j+1,i) - 2.d0 * u(k,j,i) + u(k,j-1,i)) / dyy
               dvy(k,j,i) = (v(k,j+1,i) - 2.d0 * v(k,j,i) + v(k,j-1,i)) / dyy
               dwy(k,j,i) = (w(k,j+1,i) - 2.d0 * w(k,j,i) + w(k,j-1,i)) / dyy

            elseif (j == 1) then
               duy(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j+1,i) + u(k,j+2,i)) / dyy
               dvy(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j+1,i) + v(k,j+2,i)) / dyy
               dwy(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j+1,i) + w(k,j+2,i)) / dyy

            else
               duy(k,j,i) = (u(k,j,i) - 2.d0 * u(k,j-1,i) + u(k,j-2,i)) / dyy
               dvy(k,j,i) = (v(k,j,i) - 2.d0 * v(k,j-1,i) + v(k,j-2,i)) / dyy
               dwy(k,j,i) = (w(k,j,i) - 2.d0 * w(k,j-1,i) + w(k,j-2,i)) / dyy
            endif

!           Now compute all the terms involving z derivatives, so here we
!           calculate d2u/dz2, d2v/dz2, and d2w/dz2

            if (k > 1 .and. k < nz) then
               duz(k,j,i) = (u(k+1,j,i) - 2.d0 * u(k,j,i) + u(k-1,j,i)) / dzz
               dvz(k,j,i) = (v(k+1,j,i) - 2.d0 * v(k,j,i) + v(k-1,j,i)) / dzz
               dwz(k,j,i) = (w(k+1,j,i) - 2.d0 * w(k,j,i) + w(k-1,j,i)) / dzz

            elseif (k == 1) then
               duz(k,j,i) = (u(k,j,i) - 2.d0 * u(k+1,j,i) + u(k+2,j,i)) / dzz
               dvz(k,j,i) = (v(k,j,i) - 2.d0 * v(k+1,j,i) + v(k+2,j,i)) / dzz
               dwz(k,j,i) = (w(k,j,i) - 2.d0 * w(k+1,j,i) + w(k+2,j,i)) / dzz

            else
               duz(k,j,i) = (u(k,j,i) - 2.d0 * u(k-1,j,i) + u(k-2,j,i)) / dzz
               dvz(k,j,i) = (v(k,j,i) - 2.d0 * v(k-1,j,i) + v(k-2,j,i)) / dzz
               dwz(k,j,i) = (w(k,j,i) - 2.d0 * w(k-1,j,i) + w(k-2,j,i)) / dzz
            endif

            enddo
         enddo
      enddo
      !$omp end do

   endif


!  The second block is for high-order finite difference schemes

   elseif (finite_scheme == 'high-order') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           The following is very important as it describes how we will
!           be approximating the 9 vector Laplacian terms using finite
!           differences
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
!           still use a centered difference schemes, but of lower accuracy,
!           i.e. p = 2 or p = 4. When at the boundaries of the grid,
!
!           i = 1, nx or
!           j = 1, ny or
!           k = 1, nz
!
!           use either a forward of backward difference scheme, both with
!           p = 6
!
!           First compute all the terms involving x derivatives, so here we
!           calculate d2u/dx2, d2v/dx2, and d2w/dx2

            if (i > 3 .and. i < nx - 2) then
               dux(k,j,i) = (b3 * u(k,j,i-3) + b2 * u(k,j,i-2) + &
                             b1 * u(k,j,i-1) + b0 * u(k,j,i) + &
                             b1 * u(k,j,i+1) + b2 * u(k,j,i+2) + &
                             b3 * u(k,j,i+3)) / dxx
               dvx(k,j,i) = (b3 * v(k,j,i-3) + b2 * v(k,j,i-2) + &
                             b1 * v(k,j,i-1) + b0 * v(k,j,i) + &
                             b1 * v(k,j,i+1) + b2 * v(k,j,i+2) + &
                             b3 * v(k,j,i+3)) / dxx
               dwx(k,j,i) = (b3 * w(k,j,i-3) + b2 * w(k,j,i-2) + &
                             b1 * w(k,j,i-1) + b0 * w(k,j,i) + &
                             b1 * w(k,j,i+1) + b2 * w(k,j,i+2) + &
                             b3 * w(k,j,i+3)) / dxx

            elseif (i > 2 .and. i < nx - 1) then
               dux(k,j,i) = (a2 * u(k,j,i-2) + a1 * u(k,j,i-1) + &
                             a0 * u(k,j,i) + a1 * u(k,j,i+1) + &
                             a2 * u(k,j,i+2)) / dxx
               dvx(k,j,i) = (a2 * v(k,j,i-2) + a1 * v(k,j,i-1) + &
                             a0 * v(k,j,i) + a1 * v(k,j,i+1) + &
                             a2 * v(k,j,i+2)) / dxx
               dwx(k,j,i) = (a2 * w(k,j,i-2) + a1 * w(k,j,i-1) + &
                             a0 * w(k,j,i) + a1 * w(k,j,i+1) + &
                             a2 * w(k,j,i+2)) / dxx

            elseif (i > 1 .and. i < nx) then
               dux(k,j,i) = (u(k,j,i+1) - 2.d0 * u(k,j,i) + u(k,j,i-1)) / dxx
               dvx(k,j,i) = (v(k,j,i+1) - 2.d0 * v(k,j,i) + v(k,j,i-1)) / dxx
               dwx(k,j,i) = (w(k,j,i+1) - 2.d0 * w(k,j,i) + w(k,j,i-1)) / dxx

            elseif (i == 1) then
               dux(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i+1) + &
                             c2 * u(k,j,i+2) + c3 * u(k,j,i+3) + &
                             c4 * u(k,j,i+4) + c5 * u(k,j,i+5) + &
                             c6 * u(k,j,i+6) + c7 * u(k,j,i+7)) / dxx
               dvx(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j,i+1) + &
                             c2 * v(k,j,i+2) + c3 * v(k,j,i+3) + &
                             c4 * v(k,j,i+4) + c5 * v(k,j,i+5) + &
                             c6 * v(k,j,i+6) + c7 * v(k,j,i+7)) / dxx
               dwx(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j,i+1) + &
                             c2 * w(k,j,i+2) + c3 * w(k,j,i+3) + &
                             c4 * w(k,j,i+4) + c5 * w(k,j,i+5) + &
                             c6 * w(k,j,i+6) + c7 * w(k,j,i+7)) / dxx

            else
               dux(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j,i-1) + &
                             c2 * u(k,j,i-2) + c3 * u(k,j,i-3) + &
                             c4 * u(k,j,i-4) + c5 * u(k,j,i-5) + &
                             c6 * u(k,j,i-6) + c7 * u(k,j,i-7)) / dxx
               dvx(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j,i-1) + &
                             c2 * v(k,j,i-2) + c3 * v(k,j,i-3) + &
                             c4 * v(k,j,i-4) + c5 * v(k,j,i-5) + &
                             c6 * v(k,j,i-6) + c7 * v(k,j,i-7)) / dxx
               dwx(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j,i-1) + &
                             c2 * w(k,j,i-2) + c3 * w(k,j,i-3) + &
                             c4 * w(k,j,i-4) + c5 * w(k,j,i-5) + &
                             c6 * w(k,j,i-6) + c7 * w(k,j,i-7)) / dxx
            endif

!           Now compute all the terms involving y derivatives, so here we
!           calculate d2u/dy2, d2v/dy2, and d2w/dy2

            if (j > 3 .and. j < nx - 2) then
               duy(k,j,i) = (b3 * u(k,j-3,i) + b2 * u(k,j-2,i) + &
                             b1 * u(k,j-1,i) + b0 * u(k,j,i) + &
                             b1 * u(k,j+1,i) + b2 * u(k,j+2,i) + &
                             b3 * u(k,j+3,i)) / dyy
               dvy(k,j,i) = (b3 * v(k,j-3,i) + b2 * v(k,j-2,i) + &
                             b1 * v(k,j-1,i) + b0 * v(k,j,i) + &
                             b1 * v(k,j+1,i) + b2 * v(k,j+2,i) + &
                             b3 * v(k,j+3,i)) / dyy
               dwy(k,j,i) = (b3 * w(k,j-3,i) + b2 * w(k,j-2,i) + &
                             b1 * w(k,j-1,i) + b0 * w(k,j,i) + &
                             b1 * w(k,j+1,i) + b2 * w(k,j+2,i) + &
                             b3 * w(k,j+3,i)) / dyy

            elseif (j > 2 .and. j < nx - 1) then
               duy(k,j,i) = (a2 * u(k,j-2,i) + a1 * u(k,j-1,i) + &
                             a0 * u(k,j,i) + a1 * u(k,j+1,i) + &
                             a2 * u(k,j+2,i)) / dyy
               dvy(k,j,i) = (a2 * v(k,j-2,i) + a1 * v(k,j-1,i) + &
                             a0 * v(k,j,i) + a1 * v(k,j+1,i) + &
                             a2 * v(k,j+2,i)) / dyy
               dwy(k,j,i) = (a2 * w(k,j-2,i) + a1 * w(k,j-1,i) + &
                             a0 * w(k,j,i) + a1 * w(k,j+1,i) + &
                             a2 * w(k,j+2,i)) / dyy

            elseif (j > 1 .and. j < ny) then
               duy(k,j,i) = (u(k,j+1,i) - 2.d0 * u(k,j,i) + u(k,j-1,i)) / dyy
               dvy(k,j,i) = (v(k,j+1,i) - 2.d0 * v(k,j,i) + v(k,j-1,i)) / dyy
               dwy(k,j,i) = (w(k,j+1,i) - 2.d0 * w(k,j,i) + w(k,j-1,i)) / dyy

            elseif (j == 1) then
               duy(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j+1,i) + &
                             c2 * u(k,j+2,i) + c3 * u(k,j+3,i) + &
                             c4 * u(k,j+4,i) + c5 * u(k,j+5,i) + &
                             c6 * u(k,j+6,i) + c7 * u(k,j+7,i)) / dyy
               dvy(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j+1,i) + &
                             c2 * v(k,j+2,i) + c3 * v(k,j+3,i) + &
                             c4 * v(k,j+4,i) + c5 * v(k,j+5,i) + &
                             c6 * v(k,j+6,i) + c7 * v(k,j+7,i)) / dyy
               dwy(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j+1,i) + &
                             c2 * w(k,j+2,i) + c3 * w(k,j+3,i) + &
                             c4 * w(k,j+4,i) + c5 * w(k,j+5,i) + &
                             c6 * w(k,j+6,i) + c7 * w(k,j+7,i)) / dyy

            else
               duy(k,j,i) = (c0 * u(k,j,i) + c1 * u(k,j-1,i) + &
                             c2 * u(k,j-2,i) + c3 * u(k,j-3,i) + &
                             c4 * u(k,j-4,i) + c5 * u(k,j-5,i) + &
                             c6 * u(k,j-6,i) + c7 * u(k,j-7,i)) / dyy
               dvy(k,j,i) = (c0 * v(k,j,i) + c1 * v(k,j-1,i) + &
                             c2 * v(k,j-2,i) + c3 * v(k,j-3,i) + &
                             c4 * v(k,j-4,i) + c5 * v(k,j-5,i) + &
                             c6 * v(k,j-6,i) + c7 * v(k,j-7,i)) / dyy
               dwy(k,j,i) = (c0 * w(k,j,i) + c1 * w(k,j-1,i) + &
                             c2 * w(k,j-2,i) + c3 * w(k,j-3,i) + &
                             c4 * w(k,j-4,i) + c5 * w(k,j-5,i) + &
                             c6 * w(k,j-6,i) + c7 * w(k,j-7,i)) / dyy
            endif

!           Now compute all the terms involving z derivatives, so here we
!           calculate d2u/dz2, d2v/dz2, and d2w/dz2

            if (k > 3 .and. k < nz - 2) then
               duz(k,j,i) = (b3 * u(k-3,j,i) + b2 * u(k-2,j,i) + &
                             b1 * u(k-1,j,i) + b0 * u(k,j,i) + &
                             b1 * u(k+1,j,i) + b2 * u(k+2,j,i) + &
                             b3 * u(k+3,j,i)) / dzz
               dvz(k,j,i) = (b3 * v(k-3,j,i) + b2 * v(k-2,j,i) + &
                             b1 * v(k-1,j,i) + b0 * v(k,j,i) + &
                             b1 * v(k+1,j,i) + b2 * v(k+2,j,i) + &
                             b3 * v(k+3,j,i)) / dzz
               dwz(k,j,i) = (b3 * w(k-3,j,i) + b2 * w(k-2,j,i) + &
                             b1 * w(k-1,j,i) + b0 * w(k,j,i) + &
                             b1 * w(k+1,j,i) + b2 * w(k+2,j,i) + &
                             b3 * w(k+3,j,i)) / dzz

            elseif (k > 2 .and. k < nz - 1) then
               duz(k,j,i) = (a2 * u(k-2,j,i) + a1 * u(k-1,j,i) + &
                             a0 * u(k,j,i) + a1 * u(k+1,j,i) + &
                             a2 * u(k+2,j,i)) / dzz
               dvz(k,j,i) = (a2 * v(k-2,j,i) + a1 * v(k-1,j,i) + &
                             a0 * v(k,j,i) + a1 * v(k+1,j,i) + &
                             a2 * v(k+2,j,i)) / dzz
               dwz(k,j,i) = (a2 * w(k-2,j,i) + a1 * w(k-1,j,i) + &
                             a0 * w(k,j,i) + a1 * w(k+1,j,i) + &
                             a2 * w(k+2,j,i)) / dzz

            elseif (k > 1 .and. k < nz) then
               duz(k,j,i) = (u(k+1,j,i) - 2.d0 * u(k,j,i) + u(k-1,j,i)) / dzz
               dvz(k,j,i) = (v(k+1,j,i) - 2.d0 * v(k,j,i) + v(k-1,j,i)) / dzz
               dwz(k,j,i) = (w(k+1,j,i) - 2.d0 * w(k,j,i) + w(k-1,j,i)) / dzz

            elseif (k == 1) then
               duz(k,j,i) = (c0 * u(k,j,i) + c1 * u(k+1,j,i) + &
                             c2 * u(k+2,j,i) + c3 * u(k+3,j,i) + &
                             c4 * u(k+4,j,i) + c5 * u(k+5,j,i) + &
                             c6 * u(k+6,j,i) + c7 * u(k+7,j,i)) / dzz
               dvz(k,j,i) = (c0 * v(k,j,i) + c1 * v(k+1,j,i) + &
                             c2 * v(k+2,j,i) + c3 * v(k+3,j,i) + &
                             c4 * v(k+4,j,i) + c5 * v(k+5,j,i) + &
                             c6 * v(k+6,j,i) + c7 * v(k+7,j,i)) / dzz
               dwz(k,j,i) = (c0 * w(k,j,i) + c1 * w(k+1,j,i) + &
                             c2 * w(k+2,j,i) + c3 * w(k+3,j,i) + &
                             c4 * w(k+4,j,i) + c5 * w(k+5,j,i) + &
                             c6 * w(k+6,j,i) + c7 * w(k+7,j,i)) / dzz

            else
               duz(k,j,i) = (c0 * u(k,j,i) + c1 * u(k-1,j,i) + &
                             c2 * u(k-2,j,i) + c3 * u(k-3,j,i) + &
                             c4 * u(k-4,j,i) + c5 * u(k-5,j,i) + &
                             c6 * u(k-6,j,i) + c7 * u(k-7,j,i)) / dzz
               dvz(k,j,i) = (c0 * v(k,j,i) + c1 * v(k-1,j,i) + &
                             c2 * v(k-2,j,i) + c3 * v(k-3,j,i) + &
                             c4 * v(k-4,j,i) + c5 * v(k-5,j,i) + &
                             c6 * v(k-6,j,i) + c7 * v(k-7,j,i)) / dzz
               dwz(k,j,i) = (c0 * w(k,j,i) + c1 * w(k-1,j,i) + &
                             c2 * w(k-2,j,i) + c3 * w(k-3,j,i) + &
                             c4 * w(k-4,j,i) + c5 * w(k-5,j,i) + &
                             c6 * w(k-6,j,i) + c7 * w(k-7,j,i)) / dzz
            endif

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

   return

end subroutine full_wind
