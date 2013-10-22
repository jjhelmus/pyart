subroutine biggerstaff_class(z, ze, sclass, dx, dy, dz, work_level, lapse_lim, &
                             grad_lim, z_melt, layer, win_size, fill_val, &
                             nx, ny, nz, bclass)

   implicit none

   integer(kind=4) :: nx, ny, nz
   integer(kind=4) :: win_size
   real(kind=8), intent(in) :: dx, dy, dz, work_level, lapse_lim, grad_lim
   real(kind=8), intent(in) :: z_melt, layer, fill_val
   real(kind=8), intent(in), dimension(nz) :: z
   real(kind=8), intent(in), dimension(nz,ny,nx) :: ze
   real(kind=8), intent(in), dimension(ny,nx) :: sclass
   real(kind=8), intent(out), dimension(ny,nx) :: bclass

   !!! Local variables !!!

   ! Algorithm variables
   real(kind=8), dimension(ny,nx) :: dzedx, dzedy, dzedz, grad_ze

   ! Index variables
   integer(kind=4) :: i, j, k, l, m, imin, imax, jmin, jmax, kmin, kmax
   integer(kind=4) :: kml, kze

   ! Math placeholder variables
   logical, dimension(nz,ny,nx) :: ze_m

   ! Finite difference variables
   real(kind=8), parameter :: c1=1.d0/60.d0, c2=3.d0/20.d0, c3=3.d0/4.d0
   real(kind=8), parameter :: o0=-49.d0/20.d0, o1=6.d0, o2=-15.d0/2.d0
   real(kind=8), parameter :: o3=20.d0/3.d0, o4=-15.d0/4.d0, o5=6.d0/5.d0
   real(kind=8), parameter :: o6=-1.d0/6.d0

   !!!!!!!!!!!!!!!!!!!!!!!

   !f2py intent(in, optional) nx, ny, nz
   !f2py intent(in) z, ze, sclass, dx, dy, dz, sep_alt, lapse_lim, grad_lim
   !f2py intent(in) z_melt, layer, win_size, fill_val
   !f2py intent(out) bclass

   ! Initialize Biggerstaff echo classification array
   bclass = sclass

   ! Fill appropriate variables
   dzedz = fill_val
   dzedx = fill_val
   dzedy = fill_val

   ! Create reflectivity mask which highlights invalid grid points
   ze_m = ze /= fill_val

   ! Get index of melting level
   kml = minloc(abs(z - z_melt), dim=1)

   ! Compute reflectivity lapse rate in (dB/km) within specified layer
   ! Reflectivity lapse rate is defined as the decrease in reflectivity
   ! with increasing height within the specified layer
   do i = 1,nx
      do j = 1,ny
         ! First check if column is valid for echo classification
         ! by using Steiner classification
         if (sclass(j,i) /= fill_val) then
            ! Get index of max reflectivity in column
            ! Get index corresponding to top of specified layer
            ! above max reflectivity
            kmin = maxloc(ze(:,j,i), 1, ze_m(:,j,i))
            kmax = min(nz, int(kmin + layer/dz))
            ! Compute reflectivity lapse rate within layer above
            ! maximum reflectivity in (dB/m)
            if (ze_m(kmax,j,i)) then
               dzedz(j,i) = (ze(kmin,j,i) - ze(kmax,j,i)) / layer
            endif
         endif
      enddo
   enddo
   where (dzedz /= fill_val) dzedz = 1000.d0 * dzedz

   ! Get index of working level
   k = minloc(abs(z - work_level), dim=1)

   ! Compute dzedx in (dB/m) at the working level
   do j = 1,ny
      do i = 1,nx
         if (i-1 >= 1 .and. i+1 <= nx) then
            ! 2nd-order accurate centered difference
            if (all(ze_m(k,j,i-1:i+1), dim=1)) then
               dzedx(j,i) = (ze(k,j,i+1) - ze(k,j,i-1)) / (2.d0 * dx)
            endif
         elseif (i-1 < 1) then
            ! 1st-order accurate forward difference
            if (all(ze_m(k,j,i:i+1), dim=1)) then
               dzedx(j,i) = (ze(k,j,i+1) - ze(k,j,i)) / dx
            endif
         else
            ! 1st-order accurate backward difference
            if (all(ze_m(k,j,i-1:i), dim=1)) then
               dzedx(j,i) = (ze(k,j,i) - ze(k,j,i-1)) / dx
            endif
         endif
      enddo
   enddo

   ! Compute dzedy in (dB/m) at the working level
   do i = 1,nx
      do j = 1,ny
         if (j-1 >= 1 .and. j+1 <= ny) then
            ! 2nd-order accurate centered difference
            if (all(ze_m(k,j-1:j+1,i), dim=1)) then
               dzedy(j,i) = (ze(k,j+1,i) - ze(k,j-1,i)) / (2.d0 * dy)
            endif
         elseif (j-1 < 1) then
            ! 1st-order accurate forward difference
            if (all(ze_m(k,j:j+1,i), dim=1)) then
               dzedy(j,i) = (ze(k,j+1,i) - ze(k,j,i)) / dy
            endif
         else
            ! 1st-order accurate backward difference
            if (all(ze_m(k,j-1:j,i), dim=1)) then
               dzedy(j,i) = (ze(k,j,i) - ze(k,j-1,i)) / dy
            endif
         endif
      enddo
   enddo

   ! Compute magnitude of 2-D horizontal reflectivity
   ! gradient in (dB/km)
   ! Fill gradient where we don't have valid reflectivity data
   grad_ze = 1000.d0 * sqrt(dzedx**2 + dzedy**2) ! in (dB/km)
   where (dzedx == fill_val) grad_ze = fill_val
   where (dzedy == fill_val) grad_ze = fill_val

   ! Use criteria to determine whether Steiner classifications
   ! should be changed or not

   ! Apply windowing procedure to classified field
   do i = 1,nx
      imin = max(1, i - win_size)
      imax = min(nx, i + win_size)
      do j = 1,ny
         jmin = max(1, j - win_size)
         jmax = min(ny, j + win_size)
         do l = imin,imax
            do m = jmin,jmax
            enddo
         enddo
      enddo
   enddo

end subroutine biggerstaff_class
