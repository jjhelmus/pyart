!  Module: background.f90


subroutine wind_cost(u, v, w, ub, vb, wb, wgt_ub, wgt_vb, wgt_wb, &
                     wgt_w0, fill_value, proc, nx, ny, nz, Jb)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   real(kind=8), intent(in)                       :: wgt_ub, wgt_vb, &
                                                     wgt_wb, wgt_w0, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w, &
                                                     ub, vb, wb
   real(kind=8), intent(out)                      :: Jb


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py real(kind=8), intent(in)              :: wgt_ub, wgt_vb, wgt_wb
   !f2py real(kind=8), intent(in)              :: wgt_w0, fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w, ub, vb, wb
   !f2py real(kind=8), intent(out)             :: Jb

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the background cost Jb, which is
!  given by,
!
!  Jb = 0.5 * [ sum( wgt_ub * (u - ub)**2 ) + sum( wgt_vb * (v - vb)**2 ) +
!             + sum( wgt_wb * (w - wb)**2 ) ]
!
!  where the summations are over the N Cartesian grid points. Note how the
!  background cost is a scalar value

   Jb = 0.d0

   !$omp parallel num_threads(proc)

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the value of the background cost by summing all of its
!        values at each grid point

         if (k == 1) then
            Jb = Jb + 0.5d0 * (wgt_ub * (u(k,j,i) - ub(k,j,i))**2 + &
                               wgt_vb * (v(k,j,i) - vb(k,j,i))**2 + &
                               wgt_wb * (w(k,j,i) - wb(k,j,i))**2 + &
                               wgt_w0 * w(k,j,i)**2)

         else
            Jb = Jb + 0.5d0 * (wgt_ub * (u(k,j,i) - ub(k,j,i))**2 + &
                               wgt_vb * (v(k,j,i) - vb(k,j,i))**2 + &
                               wgt_wb * (w(k,j,i) - wb(k,j,i))**2)
         endif

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_cost


subroutine wind_grad(u, v, w, ub, vb, wb, wgt_ub, wgt_vb, wgt_wb, wgt_w0, &
                     fill_value, proc, nx, ny, nz, dJbu, dJbv, dJbw)

   implicit none

   integer(kind=4), intent(in)                    :: nx, ny, nz, proc
   real(kind=8), intent(in)                       :: wgt_ub, wgt_vb, &
                                                     wgt_wb, wgt_w0, &
                                                     fill_value
   real(kind=8), intent(in), dimension(nz,ny,nx)  :: u, v, w, &
                                                     ub, vb, wb
   real(kind=8), intent(out), dimension(nz,ny,nx) :: dJbu, dJbv, dJbw


!  Define local variables ====================================================

   integer(kind=4) :: i, j, k

!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nx, ny, nz
   !f2py integer(kind=4), intent(in)           :: proc
   !f2py real(kind=8), intent(in)              :: wgt_ub, wgt_vb, wgt_wb
   !f2py real(kind=8), intent(in)              :: wgt_w0, fill_value
   !f2py real(kind=8), intent(in)              :: u, v, w, ub, vb, wb
   !f2py real(kind=8), intent(out)             :: Jb

!  ===========================================================================


!  Recall that we are attempting to minimize a function of the form,
!
!  J = J(u1,u2,...,uN,v1,v2,...,vN,w1,w2,...,wN)
!
!  which is a function of 3N variables. Note that J is typically the sum of
!  multiple different costs, including the background cost Jb, which is
!  given by,
!
!  Jb = 0.5 * [ sum( wgt_ub * (u - ub)**2 ) + sum( wgt_vb * (v - vb)**2 ) +
!             + sum( wgt_wb * (w - wb)**2 ) ]
!
!  where the summations are over the N Cartesian grid points. Note how the
!  background cost is a scalar value
!
!  We need to compute dJ/du, dJ/dv, and dJ/dw, since a minimum in J
!  corresponds with these 3 derivatives vanishing. Therefore, we need to
!  compute dJb/du, dJb/dv, and dJb/dw. Each of these terms will eventually
!  need to be vectors of length N since,
!
!  dJ/d(u,v,w) = (dJ/du1,...,dJ/duN,dJ/dv1,...,dJ/dvN,dJ/dw1,...dJ/dwN)
!
!  We minimize J in the vector space (1-D) as shown above, but we will
!  initially compute the gradient of Jb in the so-called grid space (3-D),
!  since this is usually the most natural way. These 3-D arrays will
!  eventually have to be permuted to vectors outside of this subroutine

   !$omp parallel num_threads(proc)

   !$omp do
   do i = 1, nx
      do j = 1, ny
         do k = 1, nz

!        Compute the gradient of the background cost with respect to the
!        3 control variables (u,v,w), which means we need to compute dJb/du,
!        dJb/dv, and dJb/dw. These terms are easily derived from Jb,
!
!        dJb/du = wgt_ub * (u - ub) for all N
!        dJb/dv = wgt_vb * (v - vb) for all N
!        dJb/dw = wgt_wb * (w - wb) for all N

         if (k == 1) then
            dJbu(k,j,i) = wgt_ub * (u(k,j,i) - ub(k,j,i))
            dJbv(k,j,i) = wgt_vb * (v(k,j,i) - vb(k,j,i))
            dJbw(k,j,i) = wgt_wb * (w(k,j,i) - wb(k,j,i)) + wgt_w0 * w(k,j,i)

         else
            dJbu(k,j,i) = wgt_ub * (u(k,j,i) - ub(k,j,i))
            dJbv(k,j,i) = wgt_vb * (v(k,j,i) - vb(k,j,i))
            dJbw(k,j,i) = wgt_wb * (w(k,j,i) - wb(k,j,i))
         endif

         enddo
      enddo
   enddo
   !$omp end do

   !$omp end parallel

   return

end subroutine wind_grad
