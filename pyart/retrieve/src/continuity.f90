!  Module: continuity.f90


subroutine wind_cost_orig(w, wc, wgt_c, fill_value, nx, ny, nz, Jc)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: wgt_c, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: w, wc
   real(kind=8), intent(out)                      :: Jc


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: w, wc, wgt_c, fill_value
   !f2py real(kind=8), intent(out)             :: Jc

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the anelastic air mass continuity cost
!  Jc, which in this case is given by,
!
!  Jc = 0.5 * sum( wgt_c * (w - wc)**2 )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  continuity cost is a scalar value

   Jc = 0.d0

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the value of the continuity cost by summing all of
!        its values at each grid point

         Jc = Jc + 0.5d0 * wgt_c * (w(k,j,i) - wc(k,j,i))**2

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_cost_orig


subroutine wind_grad_orig(w, wc, wgt_c, fill_value, nx, ny, nz, &
                          dJcu, dJcv, dJcw)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: wgt_c, fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: w, wc
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJcu, dJcv, dJcw


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: w, wc, wgt_c, fill_value
   !f2py real(kind=8), intent(out)             :: dJcu, dJcv, dJcw

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the anelastic air mass continuity cost
!  Jc, which in this case is given by,
!
!  Jc = 0.5 * sum( wgt_c * (w - wc)**2 )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  continuity cost is a scalar value
!
!  We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J
!  corresponds with these 3 derivatives vanishing. Therefore, we need to
!  compute dJc/du, dJc/dv, and dJc/dw. Each of these terms will eventually
!  need to be vectors of length N since,
!
!  dJ/d(u,v,w) = (dJ/du1,...,dJ/duN,dJ/dv1,...,dJ/dvN,dJ/dw1,...dJ/dwN)
!
!  We minimize J in the vector space (1-D) as shown above, but we will
!  initially compute the gradient of Jc in the so-called grid space (3-D),
!  since this is usually the most natural way. These 3-D arrays will
!  eventually have to be permuted to vectors outside of this subroutine
!
!  Now given the definition of Jc above, we quickly see that it has no
!  dependence on u or v, which greatly simplifies the problem since dJc/du
!  and dJc/dv will both be 0,
!
!  dJc/du = (0,...,0) for all (u1,u2,...,uN)
!  dJc/dv = (0,...,0) for all (v1,v2,...,vN)

  !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the gradient of the continuity cost with respect to the
!        3 control variables (u,v,w), which means we need to compute dJc/du,
!        dJc/dv, and dJc/dw. However, since the continuity cost has no
!        dependence on u or v, then we have
!
!        dJc/du = 0 for all N
!        dJc/dv = 0 for all N
!
!        Furthemore, note how dJc/dw is easily derived from Jc,
!
!        dJc/dw = wgt_c * (w - wc) for all N

         dJcu(k,j,i) = 0.d0
         dJcv(k,j,i) = 0.d0
         dJcw(k,j,i) = wgt_c * (w(k,j,i) - wc(k,j,i))

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_grad_orig


subroutine wind_cost_potvin(w, du, dv, dw, rho, drho, wgt_c, fill_value, &
                            nx, ny, nz, Jc)
                
   implicit none
   
   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: wgt_c, fill_value
   real(kind=8), intent(in), dimension(nz)        :: rho, drho
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: w, du, dv, dw
   real(kind=8), intent(out)                      :: Jc
   

!  Local variables ===========================================================

   real(kind=8)            :: D

   integer(kind=4)         :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================
   
   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: wgt_c, fill_value
   !f2py real(kind=8), intent(in)              :: w, du, dv, dw, rho, drho
   !f2py real(kind=8), intent(out)             :: Jc
   
!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the anelastic air mass continuity cost
!  Jc, which in this case is given by,
!
!  Jc = 0.5 * sum( wgt_c * [ du/dx + dv/dy + dw/dz + w/rho * drho/dz]**2 )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  continuity cost is a scalar value

   Jc = 0.d0

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the value of the continuity cost by summing all of its
!        values at each grid point
!
!        First we will calculate the main term found in Jc,
!
!        D = du/dx + dv/dy + dw/dz + w/rho * drho/dz
!
!        which then needs to be squared and summed with the rest of the
!        grid points

         D = du(k,j,i) + dv(k,j,i) + dw(k,j,i) + w(k,j,i) * drho(k) / rho(k)
                
         Jc = Jc + 0.5d0 * wgt_c * D**2

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return
                
end subroutine wind_cost_potvin


subroutine wind_grad_potvin(w, du, dv, dw, rho, drho, wgt_c, dx, dy, dz, &
                            finite_scheme, fill_value, nx, ny, nz, &
                            dJcu, dJcv, dJcw)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   character(len=16), intent(in)                  :: finite_scheme
   real(kind=8), intent(in)                       :: wgt_c, dx, dy, dz, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nz)        :: rho, drho
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: w, du, dv, dw
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJcu, dJcv, dJcw


!  Local variables ===========================================================

   real(kind=8), parameter           :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                                        b1=3.d0/4.d0, b2=3.d0/20.d0, &
                                        b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                                        c1=6.d0, c2=-15.d0/2.d0, &
                                        c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                                        c5=6.d0/5.d0, c6=-1.d0/6.d0

   real(kind=8), dimension(nz,ny,nx) :: D

   integer(kind=4)                   :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: wgt_c, dx, dy, dz
   !f2py real(kind=8), intent(in)              :: rho, drho, fill_value
   !f2py real(kind=8), intent(in)              :: w, du, dv, dw
   !f2py real(kind=8), intent(out)             :: dJcu, dJcv, dJcw

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the anelastic air mass continuity cost
!  Jc, which in this case is given by,
!
!  Jc = 0.5 * sum( wgt_c * [ du/dx + dv/dy + dw/dz + w/rho * drho/dz]**2 )
!
!  where the summation is over the N Cartesian grid points. Note how the
!  continuity cost is a scalar value.
!
!  We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J
!  corresponds with these 3 derivatives vanishing. Therefore, we need to
!  compute dJc/du, dJc/dv, and dJc/dw. Each of these terms will eventually
!  need to be vectors of length N since,
!
!  dJ/d(u,v,w) = (dJ/du1,...,dJ/duN,dJ/dv1,...,dJ/dvN,dJ/dw1,...dJ/dwN)
!
!  We minimize J in the vector space (1-D) as shown above, but we will
!  initially compute the gradient of Jc in the so-called grid space (3-D),
!  since this is usually the most natural way. These 3-D arrays will
!  eventually have to be permuted to vectors outside of this subroutine
!
!  The partial derivatives in the definition of Jc above must be
!  approximated by finite differences. This means that an analytical
!  solution to dJc/du, dJc/dv and dJc/dw does not exist. Each one of these
!  terms requires an in-depth analysis as to how it is computed. In
!  particular, it requires knowledge of the underlying finite difference
!  schemes used to approximate the partial derivatives since the
!  (u1,...,uN,v1,...,vN,w1,...,wN) terms may have influence on surrounding
!  grid points

   !$omp parallel
   
!  Calculate the main term found in Jc,
!
!  D = du/dx + dv/dy + dw/dz + w/rho * drho/dz
!
!  We need this term because Jc is proportional to the square of this term,
!  which means the gradient of Jc with respect to the 3 control variables 
!  will also depend on this term via the chain rule 

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz
         
         D(k,j,i) = du(k,j,i) + dv(k,j,i) + dw(k,j,i) + &
                    w(k,j,i) * drho(k) / rho(k)
         
         enddo
      enddo
   enddo
   !$omp end do

!  First block is for when basic finite difference schemes have been used to
!  compute the 3-D wind divergence

   if (finite_scheme == 'basic') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz
            
!           Compute the gradient of the continuity cost with respect to
!           the 3 control variables (u,v,w), which means we need to compute
!           dJc/du, dJc/dv, and dJc/dw
!
!           First compute dJc/du, which only depends on the du/dx term since
!           the other terms are independent of u

            if (i > 2 .and. i < nx - 1) then
               dJcu(k,j,i) = wgt_c * (D(k,j,i-1) - &
                                      D(k,j,i+1)) / (2.d0 * dx)

            elseif (i == 2) then
               dJcu(k,j,i) = wgt_c * (2.d0 * D(k,j,i-1) - &
                                      D(k,j,i+1)) / (2.d0 * dx)

            elseif (i == 1) then
               dJcu(k,j,i) = -wgt_c * (2.d0 * D(k,j,i) + &
                                       D(k,j,i+1)) / (2.d0 * dx)

            elseif (i == nx - 1) then
               dJcu(k,j,i) = wgt_c * (D(k,j,i-1) - 2.d0 * &
                                      D(k,j,i+1)) / (2.d0 * dx)

            else
               dJcu(k,j,i) = wgt_c * (D(k,j,i-1) + 2.d0 * &
                                      D(k,j,i)) / (2.d0 * dx)
            endif

!           Now compute dJc/dv, which only depends on the dv/dy term since
!           the other terms are independent of v

            if (j > 2 .and. j < ny - 1) then
               dJcv(k,j,i) = wgt_c * (D(k,j-1,i) - &
                                      D(k,j+1,i)) / (2.d0 * dy)

            elseif (j == 2) then
               dJcv(k,j,i) = wgt_c * (2.d0 * D(k,j-1,i) - &
                                      D(k,j+1,i)) / (2.d0 * dy)

            elseif (j == 1) then
               dJcv(k,j,i) = -wgt_c * (2.d0 * D(k,j,i) + &
                                       D(k,j+1,i)) / (2.d0 * dy)

            elseif (j == ny - 1) then
               dJcv(k,j,i) = wgt_c * (D(k,j-1,i) - 2.d0 * &
                                      D(k,j+1,i)) / (2.d0 * dy)

            else
               dJcv(k,j,i) = wgt_c * (D(k,j-1,i) + 2.d0 * &
                                      D(k,j,i)) / (2.d0 * dy)
            endif

!           Now compute dJc/dw, which only depends on the dw/dz term and
!           the w / rho * drho/dz term, since the other terms are
!           independent of w

            if (k > 2 .and. k < nz - 1) then
               dJcw(k,j,i) = wgt_c * ((D(k-1,j,i) - &
                                       D(k+1,j,i)) / (2.d0 * dz) + &
                                       D(k,j,i) * drho(k) / rho(k))

            elseif (k == 2) then
               dJcw(k,j,i) = wgt_c * ((2.d0 * D(k-1,j,i) - &
                                       D(k+1,j,i)) / (2.d0 * dz) + &
                                       D(k,j,i) * drho(k) / rho(k))

            elseif (k == 1) then
               dJcw(k,j,i) = -wgt_c * ((2.d0 * D(k,j,i) + &
                                        D(k+1,j,i)) / (2.d0 * dz) - & 
                                        D(k,j,i) * drho(k) / rho(k))

            elseif (k == nz - 1) then
               dJcw(k,j,i) = wgt_c * ((D(k-1,j,i) - 2.d0 * &
                                       D(k+1,j,i)) / (2.d0 * dz) + &
                                       D(k,j,i) * drho(k) / rho(k))
                                      

            else
               dJcw(k,j,i) = wgt_c * ((D(k-1,j,i) + 2.d0 * &
                                       D(k,j,i)) / (2.d0 * dz) + &
                                       D(k,j,i) * drho(k) / rho(k))
            endif

            enddo
         enddo
      enddo
      !$omp end do


!  The second block is for when high-order finite difference schemes have
!  been used to compute the 3-D wind divergence

   elseif (finite_scheme == 'high-order') then

      !$omp do
      do i = 1, nx
         do j = 1, ny
            do k = 1, nz

!           TODO: Add capabilites to compute gradients when high-order
!           finite difference schemes have been used

            enddo
         enddo
      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

   return

end subroutine wind_grad_potvin


subroutine integrate_up(div, rho, drhodz, dz, fill_value, nx, ny, nz, w)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: dz, fill_value
   real(kind=8), intent(in), dimension(nz)        :: rho, drhodz
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: div
   real(kind=8), intent(out), dimension(nz,ny,nx) :: w


!  Define local variables ====================================================

   real(kind=8), dimension(nz)       :: drho

   real(kind=8), dimension(nz,ny,nx) :: vel

   integer(kind=4)                   :: i, j, k

!  ===========================================================================


!  F2PY derivatives ==========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: rho, drhodz, dz
   !f2py real(kind=8), intent(in)              :: div, fill_value
   !f2py real(kind=8), intent(out)             :: w

!  ===========================================================================


!  Here we integrate the anelastic air mass continuity equation upwards
!  from the surface. The anelastic air mass continuity equation is given
!  by,
!
!  du/dx + dv/dy + dw/dz + w/rho * drho/dz = 0
!
!  One method of integrating this partial differential equation is
!  through implicit methods of finite differences

   drho = drhodz * dz

   vel = div * dz

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!           When at grid points above the bottom boundary, w is related
!           to the horizontal wind divergence, the air density and the
!           air density lapse rate

            if (k > 1) then
               w(k,j,i) = rho(k) * (w(k-1,j,i) - vel(k,j,i)) / &
                                   (rho(k) + drho(k))

!           At the bottom boundary condition, which is the surface of the
!           Earth, we require the impermeability condition, meaning that
!           w must vanish at the surface

            else
               w(k,j,i) = 0.d0

            endif

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine integrate_up


subroutine integrate_down(div, top, rho, drhodz, z, dz, fill_value, &
                          nx, ny, nz, w)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: dz, fill_value
   real(kind=8), intent(in), dimension(nz)        :: z, rho, drhodz
   real(kind=8), intent(in), dimension(ny,nx)     :: top
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: div
   real(kind=8), intent(out), dimension(nz,ny,nx) :: w

!  Local variables ===========================================================

   real(kind=8), dimension(nz)       :: drho

   real(kind=8), dimension(nz,ny,nx) :: vel

   logical, dimension(ny,nx)         :: m_top

   integer(kind=4)                   :: i, j, k


!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: z, rho, drhodz, dz
   !f2py real(kind=8), intent(in)              :: div, top, fill_value
   !f2py real(kind=8), intent(out)             :: w

!  ===========================================================================


!  Here we integrate the anelastic air mass continuity equation downwards
!  from the echo top height. The anelastic air mass continuity equation is
!  given by,
!
!  du/dx + dv/dy + dw/dz + w/rho * drho/dz = 0
!
!  One method of integrating this partial differential equation is
!  through implicit methods of finite differences

   m_top = top /= fill_value

   drho = drhodz * dz

   vel = div * dz

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny

!        If there is no valid echo top in the column, then we lack a
!        proper boundary condition to integrate the continuity equation
!        downwards from

         if (m_top(j,i)) then

            do k = nz, 1, -1

!           When at grid points below the top boundary, w is related
!           to the horizontal wind divergence, and the density and its
!           gradient

            if (z(k) < top(j,i)) then
               w(k,j,i) = rho(k) * (w(k+1,j,i) + vel(k,j,i)) / &
                                   (rho(k) - drho(k))

!           At the top boundary condition or above, which is the echo top
!           height, we require w to vanish

            else
               w(k,j,i) = 0.d0

            endif

            enddo

         else
            w(:,j,i) = 0.d0

         endif

      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine integrate_down


subroutine weight_protat(wu, wd, top, z, fill_value, nx, ny, nz, w)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz
   real(kind=8), intent(in)                       :: fill_value
   real(kind=8), intent(in), dimension(nz)        :: z
   real(kind=8), intent(in), dimension(ny,nx)     :: top
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: wu, wd
   real(kind=8), intent(out), dimension(nz,ny,nx) :: w


!  Define local variables ====================================================

   real(kind=8)                      :: f

   logical(kind=8), dimension(ny,nx) :: m_top

   integer(kind=4)                   :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py real(kind=8), intent(in)              :: z, fill_value
   !f2py real(kind=8), intent(in)              :: wu, wd, top
   !f2py real(kind=8), intent(out)             :: w

!  ===========================================================================


!  A weighted sum of the two solutions of vertical velocity estimated from
!  the upwards and downwards integration of the anelastic air mass continuity
!  equation is used to estimate the true vertical velocity in the column
!
!  The weighed sum is designed to give more weight to the solution which has
!  the lower uncertainty and error accumulation at a specific height. Closer
!  to the lower boundary condition, more weight is given to the upwards
!  integration estimate, while closer to the top boundary condition, more
!  weight is given to the downwards integration estimate

   m_top = top /= fill_value

   !$omp parallel

   !$omp do
   do i = 1, nx
      do j = 1, ny

!        If the column does not have a valid echo top, then it also does
!        not have a valid echo base, and the column is therefore poorly
!        constrained

         if (m_top(j,i)) then

            do k = 1, nz

!           At or below the echo top height, we weight the upwards and
!           downwards integrations according to the current height
!           in the column

            if (z(k) <= top(j,i)) then
               f = z(k) / top(j,i)
               w(k,j,i) = (1.d0 - f) * wu(k,j,i) + f * wd(k,j,i)

!           Above the echo top height, we expect w to vanish

            else
               w(k,j,i) = 0.d0

            endif

            enddo

         else
            w(:,j,i) = 0.d0

         endif

      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine weight_protat
