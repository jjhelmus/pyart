!  Module: gradient.f90


subroutine density1d(rho, z, finite_scheme, fill_value, nz, drho)

   implicit none

   integer(kind=4), intent(in)              :: nz
   character(len=16), intent(in)            :: finite_scheme
   real(kind=8), intent(in)                 :: fill_value
   real(kind=8), intent(in), dimension(nz)  :: rho, z
   real(kind=8), intent(out), dimension(nz) :: drho

!  Define local variables ====================================================

   real(kind=8), parameter :: a1=2.d0/3.d0, a2=1.d0/12.d0, &
                              b1=3.d0/4.d0, b2=3.d0/20.d0, &
                              b3=-1.d0/60.d0, c0=-49.d0/20.d0, &
                              c1=6.d0, c2=-15.d0/2.d0, &
                              c3=20.d0/3.d0, c4=-15.d0/4.d0, &
                              c5=6.d0/5.d0, c6=-1.d0/6.d0

   integer(kind=4)         :: k


!  ===========================================================================


!  F2PY directives ===========================================================

   !f2py integer(kind=4), optional, intent(in) :: nz
   !f2py character(len=16), intent(in)         :: finite_scheme
   !f2py real(kind=8), intent(in)              :: rho, z, fill_value
   !f2py real(kind=8), intent(out)             :: drho

!  ===========================================================================


   drho = fill_value

   !$omp parallel if(nz > 1000)

!  The first block is for the basic finite difference schemes

   if (finite_scheme == 'basic') then

      !$omp do
      do k = 1, nz

!        When in the interior of the grid,
!
!        k = [2, nz-1]
!
!        use a centered difference scheme with p = 2. When at the boundaries
!        of the grid,
!
!        k = 1, nz
!
!        use either a forward or backward difference scheme, both with p = 1

         if (k > 1 .and. k < nz) then
            drho(k) = (rho(k+1) - rho(k-1)) / (z(k+1) - z(k-1))

         elseif (k == 1) then
            drho(k) = (rho(k+1) - rho(k)) / (z(k+1) - z(k))

         else
            drho(k) = (rho(k) - rho(k-1)) / (z(k) - z(k-1))
         endif

      enddo
      !$omp end do

!  The second block is for the high-order finite difference schemes

   elseif (finite_scheme == 'high-order') then

      !$omp do
      do k = 1, nz

!        When in the interior of the grid,
!
!        k = [4, nz-2]
!
!        use a centered difference scheme with p = 6. When near the
!        boundaries of the grid,
!
!        k = 2, 3
!
!        still use centered difference schemes, but of lower accuracy, i.e.
!        p = 2 or p = 4. When at the boundaries of the grid,
!
!        k = 1, nz
!
!        use either a forward or backward difference scheme, both with p = 6

         if (k > 3 .and. k < nz - 2) then
            drho(k) = (b3 * rho(k-3) + b2 * rho(k-2) - &
                       b1 * rho(k-1) + b1 * rho(k+1) - &
                       b2 * rho(k+2) - b3 * rho(k+3))

         elseif (k > 2 .and. k < nz - 1) then
            drho(k) = (a2 * rho(k-2) - a1 * rho(k-1) + &
                       a1 * rho(k+1) - a2 * rho(k+2))

         elseif (k > 1 .and. k < nz) then
            drho(k) = (rho(k+1) - rho(k-1))

         elseif (k == 1) then
            drho(k) = (c0 * rho(k) + c1 * rho(k+1) + &
                       c2 * rho(k+2) + c3 * rho(k+3) + &
                       c4 * rho(k+4) + c5 * rho(k+5) + &
                       c6 * rho(k+6))

         else
            drho(k) = (-c0 * rho(k) - c1 * rho(k-1) - &
                        c2 * rho(k-2) - c3 * rho(k-3) - &
                        c4 * rho(k-4) - c5 * rho(k-5) - &
                        c6 * rho(k-6))
         endif

      enddo
      !$omp end do

   else

      stop

   endif

   !$omp end parallel

  return

end subroutine density1d
