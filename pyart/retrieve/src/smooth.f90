!  Module: smooth.f90


subroutine wind_cost_potvin(dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, &
                            wgt_s1, wgt_s2, wgt_s3, wgt_s4, fill_value, &
                            nx, ny, nz, Js)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: wgt_s1, wgt_s2, &
                                                     wgt_s3, wgt_s4, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: dux, duy, duz, &
                                                     dvx, dvy, dvz, &
                                                     dwx, dwy, dwz
   real(kind=8), intent(out)                      :: Js


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  =============================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: wgt_s1, wgt_s2
   !f2py real(kind=8), intent(in)              :: wgt_s3, wgt_s4, fill_value
   !f2py real(kind=8), intent(in)              :: dux, duy, duz, dvx, dvy
   !f2py real(kind=8), intent(in)              :: dvz, dwx, dwy, dwz
   !f2py real(kind=8), intent(out)             :: Js

!=============================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the smoothness cost Js, which in this
!  case is given by,
!
!  Js = 0.5 * sum( wgt_s1 * [ (d2u/dx2)**2 + (d2v/dx2)**2 +
!                           + (d2u/dy2)**2 + (d2v/dy2)**2 ] +
!                  wgt_s2 * [ (d2u/dz2)**2 + (d2v/dz2)**2 ] +
!                  wgt_s3 * [ (d2w/dx2)**2 + (d2w/dy2)**2 ] +
!                  wgt_s4 * [ (d2w/dz2)**2 ] )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  smoothness cost is a scalar value
   Js = 0.d0

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the value of the smoothness cost by summing all of its
!        values at each grid point
         Js = Js + 0.5d0 * (wgt_s1 * (dux(k,j,i)**2 + duy(k,j,i)**2 + &
                                      dvx(k,j,i)**2 + dvy(k,j,i)**2) + &
                            wgt_s2 * (duz(k,j,i)**2 + dvz(k,j,i)**2) + &
                            wgt_s3 * (dwx(k,j,i)**2 + dwy(k,j,i)**2) + &
                            wgt_s4 * (dwz(k,j,i)**2))

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_cost_potvin


subroutine wind_gradient_potvin(dux, duy, duz, dvx, dvy, dvz, dwx, dwy, dwz, &
                                wgt_s1, wgt_s2, wgt_s3, wgt_s4, dx, dy, dz, &
                                finite_scheme, fill_value, nx, ny, nz, &
                                dJsu, dJsv, dJsw)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: wgt_s1, wgt_s2, &
                                                     wgt_s3, wgt_s4, &
                                                     dx, dy, dz, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: dux, duy, duz, &
                                                     dvx, dvy, dvz, &
                                                     dwx, dwy, dwz
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJsu, dJsv, dJsw


!  Define local variables ====================================================

   real(kind=8)    :: dxx, dyy, dzz

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: wgt_s1, wgt_s2
   !f2py real(kind=8), intent(in)              :: wgt_s3, wgt_s4
   !f2py real(kind=8), intent(in)              :: dx, dy, dz, fill_value
   !f2py real(kind=8), intent(in)              :: dux, duy, duz, dvx, dvy
   !f2py real(kind=8), intent(in)              :: dvz, dwx, dwy, dwz
   !f2py real(kind=8), intent(out)             :: dJsu, dJsv, dJsw

!=============================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the smoothness cost Js, which in this
!  case is given by,
!
!  Js = 0.5 * sum( wgt_s1 * [ (d2u/dx2)**2 + (d2v/dx2)**2 +
!                           + (d2u/dy2)**2 + (d2v/dy2)**2 ] +
!                  wgt_s2 * [ (d2u/dz2)**2 + (d2v/dz2)**2 ] +
!                  wgt_s3 * [ (d2w/dx2)**2 + (d2w/dy2)**2 ] +
!                  wgt_s4 * [ (d2w/dz2)**2 ] )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  smoothness cost is a scalar value.
!
!  We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J
!  corresponds with these 3 derivatives vanishing. Therefore, we need to
!  compute dJs/du, dJs/dv, and dJs/dw. Each of these terms will eventually
!  need to be vectors of length N since,
!
!  dJ/d(u,v,w) = (dJ/du1,...,dJ/duN,dJ/dv1,...,dJ/dvN,dJ/dw1,...dJ/dwN)
!
!  We minimize J in the vector space (1-D) as shown above, but we will
!  initially compute the gradient of Js in the so-called grid space (3-D),
!  since this is usually the most natural way. These 3-D arrays will
!  eventually have to be permuted to vectors outside of this subroutine
!
!  The partial derivatives in the definition of Js above must be
!  approximated by finite differences. This means that an analytical
!  solution to dJs/du, dJs/dv and dJs/dw does not exist. Each one of these
!  terms requires an in-depth analysis as to how it is computed. In
!  particular, it requires knowledge of the underlying finite difference
!  schemes used to approximate the partial derivatives since the
!  (u1,...,uN,v1,...,vN,w1,...,wN) terms may have influence on surrounding
!  grid points
   dJsu = 0.d0
   dJsv = 0.d0
   dJsw = 0.d0
   dxx = dx**2
   dyy = dy**2
   dzz = dz**2

   !$omp parallel

!  The first block is for when basic finite difference schemes have been
!  used to approximate the partial derivatives found in Js
   if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           Compute the gradient of the smoothness cost with respect to the
!           3 control variables, which means we need to compute dJs/du,
!           dJs/dv, and dJs/dw
!
!           Recall that these 3 terms are highly dependent on the finite
!           differences used to approximate the partial derivatives, and a
!           careful analysis must be done in order to derive the gradient
!           terms, especially near the boundaries of the grid,
!
!           i = 1, 2, 3, 4, nx-3, nx-2, nx-1, nx or
!           j = 1, 2, 3, 4, ny-3, ny-2, ny-1, ny or
!           k = 1, 2, 3, 4, nz-3, nz-2, nz-1, nz
!
!           When not near the boundaries of the grid, we have a general
!           solution, but as we get closer to the boundaries the solution
!           becomes a function of the grid point in question
!
!           First compute the contribution of all the x derivative terms
!           (d2u/dx2, d2v/dx2, and d2w/dx2) to dJs/du, dJs/dv, and dJs/dw.
!           In other words, we will focus on,
!
!           Js = 0.5 * sum( wgt_s1 * [ (d2u/dx2)**2 + (d2v/dx2)**2 ] +
!                           wgt_s3 * [ (d2w/dx2)**2 ] )
            if (i > 3 .and. i < nx - 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i-1) - 2.d0 * &
                                                     dux(k,j,i) + &
                                                     dux(k,j,i+1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i-1) - 2.d0 * &
                                                     dvx(k,j,i) + &
                                                     dvx(k,j,i+1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i-1) - 2.d0 * &
                                                     dwx(k,j,i) + &
                                                     dwx(k,j,i+1)) / dxx

            elseif (i == 3) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i-2) + &
                                                     dux(k,j,i-1) - 2.d0 * &
                                                     dux(k,j,i) + &
                                                     dux(k,j,i+1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i-2) + &
                                                     dvx(k,j,i-1) - 2.d0 * &
                                                     dvx(k,j,i) + &
                                                     dvx(k,j,i+1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i-2) + &
                                                     dwx(k,j,i-1) - 2.d0 * &
                                                     dwx(k,j,i) + &
                                                     dwx(k,j,i+1)) / dxx

            elseif (i == 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i+1) - 2.d0 * &
                                                     dux(k,j,i) - 2.d0 * &
                                                     dux(k,j,i-1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i+1) - 2.d0 * &
                                                     dvx(k,j,i) - 2.d0 * &
                                                     dvx(k,j,i-1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i+1) - 2.d0 * &
                                                     dwx(k,j,i) - 2.d0 * &
                                                     dwx(k,j,i-1)) / dxx

            elseif (i == 1) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i) + &
                                                     dux(k,j,i+1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i) + &
                                                     dvx(k,j,i+1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i) + &
                                                     dwx(k,j,i+1)) / dxx

            elseif (i == nx - 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i+2) + &
                                                     dux(k,j,i+1) - 2.d0 * &
                                                     dux(k,j,i) + &
                                                     dux(k,j,i-1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i+2) + &
                                                     dvx(k,j,i+1) - 2.d0 * &
                                                     dvx(k,j,i) + &
                                                     dvx(k,j,i-1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i+2) + &
                                                     dwx(k,j,i+1) - 2.d0 * &
                                                     dwx(k,j,i) + &
                                                     dwx(k,j,i-1)) / dxx

            elseif (i == nx - 1) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i-1) - 2.d0 * &
                                                     dux(k,j,i) - 2.d0 * &
                                                     dux(k,j,i+1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i-1) - 2.d0 * &
                                                     dvx(k,j,i) - 2.d0 * &
                                                     dvx(k,j,i+1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i-1) - 2.d0 * &
                                                     dwx(k,j,i) - 2.d0 * &
                                                     dwx(k,j,i+1)) / dxx

            else
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (dux(k,j,i) + &
                                                     dux(k,j,i-1)) / dxx
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvx(k,j,i) + &
                                                     dvx(k,j,i-1)) / dxx
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwx(k,j,i) + &
                                                     dwx(k,j,i-1)) / dxx
            endif

!           Now compute the contribution of all the y derivative terms
!           (d2u/dy2, d2v/dy2, and d2w/dy2) to dJs/du, dJs/dv, and dJs/dw.
!           In other words, we will focus on,
!
!           Js = 0.5 * sum( wgt_s1 * [ (d2u/dy2)**2 + (d2v/dy2)**2 ] +
!                           wgt_s3 * [ (d2w/dy2)**2 ] )
            if (j > 3 .and. j < ny - 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j-1,i) - 2.d0 * &
                                                     duy(k,j,i) + &
                                                     duy(k,j+1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j-1,i) - 2.d0 * &
                                                     dvy(k,j,i) + &
                                                     dvy(k,j+1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j-1,i) - 2.d0 * &
                                                     dwy(k,j,i) + &
                                                     dwy(k,j+1,i)) / dyy

            elseif (j == 3) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j-2,i) + &
                                                     duy(k,j-1,i) - 2.d0 * &
                                                     duy(k,j,i) + &
                                                     duy(k,j+1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j-2,i) + &
                                                     dvy(k,j-1,i) - 2.d0 * &
                                                     dvy(k,j,i) + &
                                                     dvy(k,j+1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j-2,i) + &
                                                     dwy(k,j-1,i) - 2.d0 * &
                                                     dwy(k,j,i) + &
                                                     dwy(k,j+1,i)) / dyy

            elseif (j == 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j+1,i) - 2.d0 * &
                                                     duy(k,j,i) - 2.d0 * &
                                                     duy(k,j-1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j+1,i) - 2.d0 * &
                                                     dvy(k,j,i) - 2.d0 * &
                                                     dvy(k,j-1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j+1,i) - 2.d0 * &
                                                     dwy(k,j,i) - 2.d0 * &
                                                     dwy(k,j-1,i)) / dyy

            elseif (j == 1) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j,i) + &
                                                     duy(k,j+1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j,i) + &
                                                     dvy(k,j+1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j,i) + &
                                                     dwy(k,j+1,i)) / dyy

            elseif (j == ny - 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j+2,i) + &
                                                     duy(k,j+1,i) - 2.d0 * &
                                                     duy(k,j,i) + &
                                                     duy(k,j-1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j+2,i) + &
                                                     dvy(k,j+1,i) - 2.d0 * &
                                                     dvy(k,j,i) + &
                                                     dvy(k,j-1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j+2,i) + &
                                                     dwy(k,j+1,i) - 2.d0 * &
                                                     dwy(k,j,i) + &
                                                     dwy(k,j-1,i)) / dyy

            elseif (j == ny - 1) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j-1,i) - 2.d0 * &
                                                     duy(k,j,i) - 2.d0 * &
                                                     duy(k,j+1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j-1,i) - 2.d0 * &
                                                     dvy(k,j,i) - 2.d0 * &
                                                     dvy(k,j+1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j-1,i) - 2.d0 * &
                                                     dwy(k,j,i) - 2.d0 * &
                                                     dwy(k,j+1,i)) / dyy

            else
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s1 * (duy(k,j,i) + &
                                                     duy(k,j-1,i)) / dyy
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s1 * (dvy(k,j,i) + &
                                                     dvy(k,j-1,i)) / dyy
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s3 * (dwy(k,j,i) + &
                                                     dwy(k,j-1,i)) / dyy
            endif

!           Now compute the contribution of all the z derivative terms
!           (d2u/dz2, d2v/dz2, and d2w/dz2) to dJs/du, dJs/dv, and dJs/dw.
!           In other words, we will focus on,
!
!           Js = 0.5 * sum( wgt_s2 * [ (d2u/dz2)**2 + (d2v/dz2)**2 ] +
!                           wgt_s4 * [ (d2w/dz2)**2 ] )
            if (k > 3 .and. k < nz - 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k-1,j,i) - 2.d0 * &
                                                     duz(k,j,i) + &
                                                     duz(k+1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k-1,j,i) - 2.d0 * &
                                                     dvz(k,j,i) + &
                                                     dvz(k+1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k-1,j,i) - 2.d0 * &
                                                     dwz(k,j,i) + &
                                                     dwz(k+1,j,i)) / dzz

            elseif (k == 3) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k-2,j,i) + &
                                                     duz(k-1,j,i) - 2.d0 * &
                                                     duz(k,j,i) + &
                                                     duz(k+1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k-2,j,i) + &
                                                     dvz(k-1,j,i) - 2.d0 * &
                                                     dvz(k,j,i) + &
                                                     dvz(k+1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k-2,j,i) + &
                                                     dwz(k-1,j,i) - 2.d0 * &
                                                     dwz(k,j,i) + &
                                                     dwz(k+1,j,i)) / dzz

            elseif (k == 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k+1,j,i) - 2.d0 * &
                                                     duz(k,j,i) - 2.d0 * &
                                                     duz(k-1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k+1,j,i) - 2.d0 * &
                                                     dvz(k,j,i) - 2.d0 * &
                                                     dvz(k-1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k+1,j,i) - 2.d0 * &
                                                     dwz(k,j,i) - 2.d0 * &
                                                     dwz(k-1,j,i)) / dzz

            elseif (k == 1) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k,j,i) + &
                                                     duz(k+1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k,j,i) + &
                                                     dvz(k+1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k,j,i) + &
                                                     dwz(k+1,j,i)) / dzz

            elseif (k == nz - 2) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k+2,j,i) + &
                                                     duz(k+1,j,i) - 2.d0 * &
                                                     duz(k,j,i) + &
                                                     duz(k-1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k+2,j,i) + &
                                                     dvz(k+1,j,i) - 2.d0 * &
                                                     dvz(k,j,i) + &
                                                     dvz(k-1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k+2,j,i) + &
                                                     dwz(k+1,j,i) - 2.d0 * &
                                                     dwz(k,j,i) + &
                                                     dwz(k-1,j,i)) / dzz

            elseif (k == nz - 1) then
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k-1,j,i) - 2.d0 * &
                                                     duz(k,j,i) - 2.d0 * &
                                                     duz(k+1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k-1,j,i) - 2.d0 * &
                                                     dvz(k,j,i) - 2.d0 * &
                                                     dvz(k+1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k-1,j,i) - 2.d0 * &
                                                     dwz(k,j,i) - 2.d0 * &
                                                     dwz(k+1,j,i)) / dzz

            else
               dJsu(k,j,i) = dJsu(k,j,i) + wgt_s2 * (duz(k,j,i) + &
                                                     duz(k-1,j,i)) / dzz
               dJsv(k,j,i) = dJsv(k,j,i) + wgt_s2 * (dvz(k,j,i) + &
                                                     dvz(k-1,j,i)) / dzz
               dJsw(k,j,i) = dJsw(k,j,i) + wgt_s4 * (dwz(k,j,i) + &
                                                     dwz(k-1,j,i)) / dzz
            endif

            enddo
         enddo
      enddo
      !$omp end do


!  The second block is for when high-order finite difference schemes have
!  been used to approximate the partial derivatives found in Js
   elseif (finite_scheme == 'high-order') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           Compute the gradient of the smoothness cost with respect to the
!           3 control variables, which means we need to compute dJs/du,
!           dJs/dv, and dJs/dw
!
!           TODO: Add capabilities to compute the gradient of the smoothness
!           cost when high-order finite difference schemes have been used
!           to approximate the partial derivatives

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

   return

end subroutine wind_gradient_potvin


subroutine wind_cost_collis(dux, duy, duz, duxy, duxz, duyz, dvx, &
                            dvy, dvz, dvxy, dvxz, dvyz, wgt_s, &
                            fill_value, nx, ny, nz, Js)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: wgt_s, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: dux, duy, duz, duxy, &
                                                     duxz, duyz, dvx, dvy, &
                                                     dvz, dvxy, dvxz, dvyz
   real(kind=8), intent(out)                      :: Js


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: wgt_s, fill_value
   !f2py real(kind=8), intent(in)              :: dux, duy, duz, duxy, duxz
   !f2py real(kind=8), intent(in)              :: duyz, dvx, dvy, dvz, dvxy
   !f2py real(kind=8), intent(in)              :: dvxz, dvyz
   !f2py real(kind=8), intent(out)             :: Js

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the smoothness cost Js, which in this
!  case is given by,
!
!  Js = 0.5 * sum( wgt_s * [ (d2u/dx2)**2 + (d2v/dx2)**2 + (d2u/dy2)**2 +
!                          + (d2v/dy2)**2 + (d2u/dz2)**2 + (d2v/dz2)**2 +
!                          + (d2u/dxdy)**2 + (d2v/dxdy)**2 + (d2u/dxdz)**2 +
!                          + (d2v/dxdz)**2 + (d2u/dydz)**2 +
!                          + (d2v/dydz)**2 ] )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  smoothness cost is a scalar value
   Js = 0.d0

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the value of the smoothness cost by summing all of its
!        values at each grid point
         Js = Js + 0.5d0 * wgt_s * (dux(k,j,i)**2 + duy(k,j,i)**2 + &
                                    duz(k,j,i)**2 + dvx(k,j,i)**2 + &
                                    dvy(k,j,i)**2 + dvz(k,j,i)**2 + &
                                    duxy(k,j,i)**2 + duxz(k,j,i)**2 + &
                                    duyz(k,j,i)**2 + dvxy(k,j,i)**2 + &
                                    dvxz(k,j,i)**2 + dvyz(k,j,i)**2)

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_cost_collis


subroutine wind_gradient_collis(dux, duy, duz, duxy, duxz, duyz, dvx, &
                                dvy, dvz, dvxy, dvxz, dvyz, wgt_s, &
                                dx, dy, dz, finite_scheme, fill_value, &
                                nx, ny, nz, dJsu, dJsv, dJsw)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: wgt_s, dx, dy, dz, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: dux, duy, duz, duxy, &
                                                     duxz, duyz, dvx, dvy, &
                                                     dvz, dvxy, dvxz, dvyz
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJsu, dJsv, dJsw


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: wgt_s, fill_value
   !f2py real(kind=8), intent(in)              :: dux, duy, duz, duxy, duxz
   !f2py real(kind=8), intent(in)              :: duyz, dvx, dvy, dvz, dvxy
   !f2py real(kind=8), intent(in)              :: dvxz, dvyz, dx, dy, dz
   !f2py real(kind=8), intent(out)             :: dJsu, dJsv, dJsw

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the smoothness cost Js, which in this
!  case is given by,
!
!  Js = 0.5 * sum( wgt_s * [ (d2u/dx2)**2 + (d2v/dx2)**2 + (d2u/dy2)**2 +
!                          + (d2v/dy2)**2 + (d2u/dz2)**2 + (d2v/dz2)**2 +
!                          + (d2u/dxdy)**2 + (d2v/dxdy)**2 + (d2u/dxdz)**2 +
!                          + (d2v/dxdz)**2 + (d2u/dydz)**2 +
!                          + (d2v/dydz)**2 ] )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  smoothness cost is a scalar value.
!
!  We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J
!  corresponds with these 3 derivatives vanishing. Therefore, we need to
!  compute dJs/du, dJs/dv, and dJs/dw. Each of these terms will eventually
!  need to be vectors of length N since,
!
!  dJ/d(u,v,w) = (dJ/du1,...,dJ/duN,dJ/dv1,...,dJ/dvN,dJ/dw1,...dJ/dwN)
!
!  We minimize J in the vector space (1-D) as shown above, but we will
!  initially compute the gradient of Js in the so-called grid space (3-D),
!  since this is usually the most natural way. These 3-D arrays will
!  eventually have to be permuted to vectors outside of this subroutine
!
!  The partial derivatives in the definition of Js above must be
!  approximated by finite differences. This means that an analytical
!  solution to dJs/du, dJs/dv and dJs/dw does not exist. Each one of these
!  terms requires an in-depth analysis as to how it is computed. In
!  particular, it requires knowledge of the underlying finite difference
!  schemes used to approximate the partial derivatives since the
!  (u1,...,uN,v1,...,vN,w1,...,wN) terms may have influence on surrounding
!  grid points
   dJsu = 0.d0
   dJsv = 0.d0
   dJsw = 0.d0

   !$omp parallel

!  The first block is for when basic finite difference schemes have been
!  used to approximate the partial derivatives found in Js

  if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           TODO: Add capabilities to compute the gradient of the smoothness
!           cost when basic finite differences have been used to approximate
!           the partial derivatives

            enddo
         enddo
      enddo
      !$omp end do


!  The second block is for when high-order finite difference schemes have
!  been used to approximate the partial derivatives found in Js

   elseif (finite_scheme == 'high-order') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
              
!           TODO: Add capabilities to compute the gradient of the smoothness
!           cost when high-order finite differences have been used to
!           approximate the partial derivatives

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

   return

end subroutine wind_gradient_collis