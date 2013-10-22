! Module: echo_class.f90

subroutine convective_radius(ze, area_relation, conv_rad)

   implicit none

   real(kind=8), intent(in) :: ze
   character(len=16), intent(in) :: area_relation
   real(kind=8), intent(out) :: conv_rad

   !!! F2PY directives !!!

   !f2py intent(in) ze, area_relation
   !f2py intent(out) conv_rad

   !!!!!!!!!!!!!!!!!!!!!!!

   if (area_relation == 'small') then
      if (ze < 30.d0) conv_rad = 1000.d0
      if (ze >= 30.d0 .and. ze < 35.d0) conv_rad = 2000.d0
      if (ze >= 35.d0 .and. ze < 40.d0) conv_rad = 3000.d0
      if (ze >= 40.d0 .and. ze < 45.d0) conv_rad = 4000.d0
      if (ze >= 45.d0) conv_rad = 5000.d0
   elseif (area_relation == 'medium') then
      if (ze < 25.d0) conv_rad = 1000.d0
      if (ze >= 25.d0 .and. ze < 30.d0) conv_rad = 2000.d0
      if (ze >= 30.d0 .and. ze < 35.d0) conv_rad = 3000.d0
      if (ze >= 35.d0 .and. ze < 40.d0) conv_rad = 4000.d0
      if (ze >= 40.d0) conv_rad = 5000.d0
   elseif (area_relation == 'large') then
      if (ze < 20.d0) conv_rad = 1000.d0
      if (ze >= 20.d0 .and. ze < 25.d0) conv_rad = 2000.d0
      if (ze >= 25.d0 .and. ze < 30.d0) conv_rad = 3000.d0
      if (ze >= 30.d0 .and. ze < 35.d0) conv_rad = 4000.d0
      if (ze >= 35.d0) conv_rad = 5000.d0
   elseif (area_relation == 'sgp') then
      if (ze < 40.d0) conv_rad = 0.d0
      if (ze >= 40.d0 .and. ze < 45.d0) conv_rad = 1000.d0
      if (ze >= 45.d0 .and. ze < 50.d0) conv_rad = 2000.d0
      if (ze >= 50.d0 .and. ze < 55.d0) conv_rad = 6000.d0 
      if (ze >= 55.d0) conv_rad = 8000.d0
   endif

 end subroutine convective_radius

 subroutine peakedness(ze, peak_relation, peak)

   implicit none

   real(kind=8), intent(in) :: ze
   character(len=16), intent(in) :: peak_relation
   real(kind=8), intent(out) :: peak

   !!! F2PY directives !!!

   !f2py intent(in) ze, peak_relation
   !f2py intent(out) peak

   !!!!!!!!!!!!!!!!!!!!!!!

   if (peak_relation == 'default') then
      if (ze < 0.d0) peak = 10.d0
      if (ze >= 0.d0 .and. ze < 42.43d0) peak = 10.d0 - ze**2 / 180.d0
      if (ze >= 42.43d0) peak = 0.d0
   elseif (peak_relation == 'sgp') then
      if (ze < 0.d0) peak = 14.d0
      if (ze >= 0.d0 .and. ze < 42.43) peak = 14.d0 - ze**2 / 180.d0
      if (ze >= 42.43) peak = 4.d0
   endif

end subroutine peakedness

subroutine steiner(x, y, z, ze, dx, dy, intense, bkg_rad, work_lev, &
                   area_relation, peak_relation, use_intense, fill_val, &
                   nx, ny, nz, sclass)

   use omp_lib

   implicit none

   integer(kind=4) :: nx, ny, nz
   logical, intent(in) :: use_intense
   real(kind=8), intent(in) :: dx, dy, intense, bkg_rad, work_lev, fill_val
   character(len=16), intent(in) :: area_relation, peak_relation
   real(kind=8), intent(in), dimension(nx) :: x
   real(kind=8), intent(in), dimension(ny) :: y
   real(kind=8), intent(in), dimension(nz) :: z
   real(kind=8), intent(in), dimension(nz,ny,nx) :: ze
   real(kind=8), intent(out), dimension(ny,nx) :: sclass

   !!! Local variables !!!

   ! Algorithm variables
   real(kind=8) :: conv_rad, peak
   real(kind=8), dimension(ny,nx) :: ze_s

   ! Index variables
   integer(kind=4) :: i, j, k, l, m, imin, imax, jmin, jmax

   ! Placeholder variables
   real(kind=8) :: n, sum_ze, bkg_ze, rad

   ! Mask variables
   logical, dimension(ny,nx) :: m_ze

   !!!!!!!!!!!!!!!!!!!!!!!

   !!! F2PY directives !!!

   !f2py intent(in, optional) nx, ny, nz
   !f2py intent(in) x, y, z, ze, col_type, dx, dy, intense, bkg_rad
   !f2py intent(in) use_intense, work_lev, area_relation
   !f2py intent(in) peak_relation, fill_val
   !f2py intent(out) sclass

   !!!!!!!!!!!!!!!!!!!!!!!

   ! Fill echo class
   sclass = fill_val

   ! Get working level index and the reflectivity
   ! slice at the working level
   k = minloc(abs(z - work_lev), dim=1)
   ze_s = ze(k,:,:)

   ! Create reflectivity mask at the working level
   m_ze = ze_s /= fill_val

   !$omp parallel

   ! Steiner algorithm
   !$omp do
   do i = 1, nx
      do j = 1, ny
         ! First we check if the grid point has not been classified, and
         ! whether it has a valid reflectivity
         if (sclass(j,i) == fill_val .and. m_ze(j,i)) then
            imin = max(1, int(i - bkg_rad / dx))
            imax = min(nx, int(i + bkg_rad / dx))
            jmin = max(1, int(j - bkg_rad / dy))
            jmax = min(ny, int(j + bkg_rad / dy))
            n = 0.d0; sum_ze = 0.d0
            do l = imin, imax
               do m = jmin, jmax
                  rad = sqrt((x(l)-x(i))**2 + (y(m)-y(j))**2) ! in (m)
                  if (rad <= bkg_rad .and. m_ze(m,l)) then
                     n = n + 1.d0
                     sum_ze = sum_ze + 10.d0**(ze_s(m,l) / 10.d0) ! in (mm^6/m^3)
                  endif
               enddo
            enddo
            bkg_ze = 10.d0 * log10(sum_ze / n) ! Background reflectivity in (dBZ)
            call convective_radius(bkg_ze, area_relation, conv_rad)
            imin = max(1, int(i - conv_rad / dx))
            imax = min(nx, int(i + conv_rad / dx))
            jmin = max(1, int(j - conv_rad / dy))
            jmax = min(ny, int(j + conv_rad / dy))
            ! Check the intensity criteria. If intensity criteria is
            ! met, surrounding grid points within the convective radius
            ! are also convective
            if (use_intense .and. ze_s(j,i) >= intense) then
               sclass(j,i) = 1.d0
               do l = imin, imax
                  do m = jmin, jmax
                     rad = sqrt((x(l)-x(i))**2 + (y(m)-y(j))**2) ! in (m)
                     if (rad <= conv_rad .and. m_ze(m,l)) then
                        sclass(m,l) = 1.d0
                     endif
                  enddo
               enddo
            ! Check the peakedness criteria. If peakedness is satisfied,
            ! surrounding grid points within the convective radius are
            ! also convective
            else
               call peakedness(bkg_ze, peak_relation, peak)
               if (ze_s(j,i) - bkg_ze >= peak) then
                  sclass(j,i) = 1.d0
                  do l = imin, imax
                     do m = jmin, jmax
                        rad = sqrt((x(l)-x(i))**2 + (y(m)-y(j))**2) ! in (m)
                        if (rad <= conv_rad .and. m_ze(m,l)) then
                           sclass(m,l) = 1.d0
                        endif
                     enddo
                  enddo
               else ! Grid point must be stratiform
                  sclass(j,i) = -1.d0
               endif
            endif
         endif
      enddo
   enddo
   !$omp end do

   !$omp end parallel

end subroutine steiner
