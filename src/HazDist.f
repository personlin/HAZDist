
      subroutine CalcDist_f (sourceType, pscorflag, nFltGrid, n1AS, iLocX, iLocY, n2AS,
     1             iFltWidth, iFlt, iMag, ystep, grid_top, RupWidth, RupLen, r_horiz, seismoDepth,
     2             fltgrid_x, fltgrid_y, fltgrid_z, fltgrid_x1, fltgrid_y1, fltgrid_z1,
     3             fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3, fltgrid_z3,
     4             fltgrid_x4, fltgrid_y4, fltgrid_z4, fltgrid_Rrup, fltgrid_Rjb, dip, dipS7,
     5             distS7, HWFlag, n1, n2, icellRupstrike, icellRupdip, hypoDepth, distJB,
     6             distRup, ZTOR, distSeismo, distepi, disthypo, dipavgd, Rx, Ry, Ry0)

c     This subroutine computes the distance parameters for the GMPEs
c     (e.g. Rjb, Rrup, ZTOR, Rx, Ry, Ry0)

      implicit none
      include 'pfrisk.h'

c     declarations passed in
      integer sourceType, pscorflag, nfltGrid(2), n1AS(MAXFLT_AS), iLocX,
     1        iLocY, n2AS(MAXFLT_AS), iFltWidth, iFlt, iMag
      real ystep, grid_top(MAX_FLT,MAX_GRID), RupWidth, RupLen, r_horiz, seismoDepth,
     1     fltgrid_x(MAXFLT_DD,MAXFLT_AS), fltgrid_y(MAXFLT_DD,MAXFLT_AS),
     2     fltgrid_z(MAXFLT_DD,MAXFLT_AS), fltgrid_x1(MAXFLT_DD,MAXFLT_AS),
     3     fltgrid_y1(MAXFLT_DD,MAXFLT_AS), fltgrid_z1(MAXFLT_DD,MAXFLT_AS),
     4     fltgrid_x2(MAXFLT_DD,MAXFLT_AS), fltgrid_y2(MAXFLT_DD,MAXFLT_AS),
     5     fltgrid_z2(MAXFLT_DD,MAXFLT_AS), fltgrid_x3(MAXFLT_DD,MAXFLT_AS),
     6     fltgrid_y3(MAXFLT_DD,MAXFLT_AS), fltgrid_z3(MAXFLT_DD,MAXFLT_AS),
     6     fltgrid_x4(MAXFLT_DD,MAXFLT_AS), fltgrid_y4(MAXFLT_DD,MAXFLT_AS),
     7     fltgrid_z4(MAXFLT_DD,MAXFLT_AS), fltgrid_Rrup(MAXFLT_DD,MAXFLT_AS),
     8     fltgrid_Rjb(MAXFLT_DD,MAXFLT_AS), dip(MAX_FLT,MAX_WIDTH,MAX_SEG),
     9     distS7(MAX_FLT,MAX_S7), dipS7(MAX_FLT,MAX_S7)

c     declarations passed out
      integer HWFlag, n1, n2, icellRupstrike, icellRupdip
      real hypoDepth, distJB, distRup, ZTOR, distSeismo, distepi, disthypo,
     1     dipavgd, Rx, Ry, Ry0

c     declarations only used within subroutine
      integer i, j, celly, cellx
      real top_grid, hypo1, rupLength2, r_horiz1, testy, testx, hypox, hypoy

c     Initialize closest distance
      distRup = 1.0e30
      distJB = 1.0e30
      distSeismo = 1.0e30
      disthypo = 1.0e30
      distepi = 1.0e30

c     Set top_grid for sourceTypes 2/3 and 4 (different)
        if (sourceType .eq. 2 .or. sourceType .eq. 3) then
          top_grid = grid_top(iFlt,1)
        else if (sourceType .eq. 4) then
          top_grid = grid_top(iFlt,iLocX)
        endif

c     Compute distances for areal sources
        if ( sourceType .eq. 2 .or. sourceType .eq. 3 .or. sourceType .eq. 4 ) then

c         Set hypoDepth for areal sources
          hypoDepth = (iLocY-0.5)*ystep + top_grid

c         Approximate correction for using point source with extended source model
C         added limit on only shallow areal sources (Hypo>30.0km)
          if ( psCorFlag .eq. 1  .and. hypoDepth .le. 30.0) then
            hypo1=hypoDepth-rupWidth/2.
          if (hypo1 .lt. 0.0) hypo1=0.0
          rupLength2 = rupLen/2.
          if ( r_horiz .gt. rupLength2 ) then
            r_horiz1 = r_horiz - 0.625*rupLength2*(1.-1./(1.5*(r_horiz/rupLength2)**(1.5)+1))
          else
            r_horiz1 = r_horiz/1.57
          endif
          else
            r_horiz1 = r_horiz
            hypo1 = hypoDepth
          endif
          distJB = r_horiz1
          distRup = sqrt( r_horiz1**2 + hypo1**2 )
          ZTOR = hypo1
          if ( hypo1 .lt. seismoDepth) then
            distSeismo = sqrt(distJB**2 + seismoDepth**2)
          else
            distSeismo = distRup
          endif

c         Set other distance values equal to corresponding distance measures.
          distepi = distJB
          disthypo = distRup
          Rx = distJB
          Ry = 0.0
          Ry0 = 0.0
C         Set HWflag = 0 for areal or grid sources.
          HWFlag = 0
          dipavgd = dip(iFlt,ifltWidth,1)
          return

c     Compute distances for SourceType 7
        elseif (sourceType .eq. 7) then
          distrup = distS7(iFlt,iMag)
          distJB  = distS7(iFlt,iMag)
          distseismo = distS7(iFlt,iMag)
          dipavgd = dipS7(iFlt,iMag)
          disthypo = distS7(iFlt,iMag)
c         Fixed Parameters
          ZTOR = 0.0
          HWFlag = 0
          hypoDepth = 8.0
          Rx = distrup

c     Otherwise use fault sources (i.e., Sourcetype = 1, 5, or 6)
        else

c         n1 is the last cell that makes up the rupture plane down dip
c         n2 is the last cell that makes up the rupture plane along strike
          n1 = nfltgrid(1)-n1AS(iLocX)+iLocY
          n2 = n2AS(iLocX)

c        Calculate hypoDepth (assumes hypocenter is in center of rupture plane)
         if (sourceType .eq. 1 .or. sourceType .eq. 5 .or. sourceType .eq. 6) then
           testy = ((n1-iLocY)+1.)/2.
           celly = int(iLocY+testy)
           testx = ((n2-iLocX)+1.)/2.
           cellx = int(iLocX+testx)
           if (int(testy) .eq. testy .and. int(testx) .eq. testx) then
             hypoDepth = fltgrid_z1(celly,cellx)
             hypox = fltgrid_x1(celly,cellx)
             hypoy = fltgrid_y1(celly,cellx)
           else if (int(testy) .eq. testy .and. int(testx) .ne. testx) then
             hypoDepth = (fltgrid_z2(celly,cellx) + fltgrid_z1(celly,cellx))/2.
             hypox = (fltgrid_x2(celly,cellx) + fltgrid_x1(celly,cellx))/2.
             hypoy = (fltgrid_y2(celly,cellx) + fltgrid_y1(celly,cellx))/2.
           else if (int(testy) .ne. testy .and. int(testx) .eq. testx) then
             hypoDepth = (fltgrid_z4(celly,cellx) + fltgrid_z1(celly,cellx))/2.
             hypox = (fltgrid_x4(celly,cellx) + fltgrid_x1(celly,cellx))/2.
             hypoy = (fltgrid_y4(celly,cellx) + fltgrid_y1(celly,cellx))/2.
           else
             hypoDepth = fltgrid_z(celly,cellx)
             hypox = fltgrid_x(celly,cellx)
             hypoy = fltgrid_y(celly,cellx)
           endif
         endif

c         Calculate Rx, Ry, Ry0, dipavgd, and HWFlag
c         Global Coordinate System 2 method
          call GC2 (iLocX, iLocY, n2, n1, fltgrid_x1, fltgrid_y1, fltgrid_z1,
     1              fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3,
     2              fltgrid_z3, fltgrid_x4, fltgrid_y4, fltgrid_z4, Rx, Ry, Ry0,
     3              HWFlag, dipavgd)

c         Calculate ZTOR
          ZTOR = fltgrid_z1(iLocY,iLocX)

c         Calculate Rrup and Rjb
C         Keep track of the fault grid cell for the closest rupture distance
C         for this rupture area since it is needed for the NGA directivity models.

          do j=iLocX,n2
            do i=iLocY,n1
              if (distRup .gt. fltgrid_Rrup(i,j)) then
                distRup = fltgrid_rRup(i,j)
                icellRupstrike = j
                icellRupdip = i
              endif

              if (distJB .gt. fltgrid_Rjb(i,j)) then
                distJB = fltgrid_Rjb(i,j)
              endif

C         Compute the DistSeismo value based on min depth being equal to seismoDepth parameter.
              if (fltgrid_z(i,j) .gt. seismoDepth) then
                if (distSeismo .gt. fltgrid_Rrup(i,j)) then
                  distSeismo = fltgrid_Rrup(i,j)
                endif
              endif
            enddo
          enddo

C         Set the Epi and Hypo distances for faults.
          distepi = sqrt(hypox**2.+hypoy**2.)
          disthypo = sqrt(hypox**2.+hypoy**2.+hypoDepth**2.)

        endif

      return
      end

c ----------------------------------------------------------------------

      Subroutine Set_MinDist (sourceType, iFlt, iFltWidth, distRup, distJB, distSeismo,
     1                        SourceDist, MinRrup_temp, MinRjb_temp, MinSeismo_temp)

      implicit none
      include 'pfrisk.h'

      integer sourceType, iFlt, iFltWidth
      real distRup, distJB, distSeismo, SourceDist(MAX_FLT,MAX_WIDTH,3),
     1     MinRrup_temp, MinRjb_temp, MinSeismo_temp

       if ( distRup .lt. SourceDist(iFlt,iFltWidth,1) ) then
         SourceDist(iFlt,iFltWidth,1)=distRup
       endif
       if ( distJB .lt. SourceDist(iFlt,iFltWidth,2) ) then
         SourceDist(iFlt,iFltWidth,2)=distJB
       endif
       if ( distSeismo .lt. SourceDist(iFlt,iFltWidth,3) ) then
         SourceDist(iFlt,iFltWidth,3)=distSeismo
       endif
       MinRrup_temp = SourceDist(iFlt,iFltWidth,1)
       MinRjb_temp = SourceDist(iFlt,iFltWidth,2)
       MinSeismo_temp = SourceDist(iFlt,iFltWidth,3)

      return
      end

c ----------------------------------------------------------------------

      Subroutine SetnRupLoc ( n1, n2, nHypoX, pHypoX, nHypoXStep,
     1                        nHypoZ, pHypoZ, nHypoZstep )

      implicit none

      integer n1, n2, nHypoX, nHypoXstep, nHypoZ, nHypoZstep
      real pHypoX, pHypoZ

C     First set up the number of hypocenter locations for a given fault rupture area
C     If there are less than 10 cells in either along strike or along dip direction
C     just use each cell. Otherwise take 10 locations along strike and dip

      if (n2 .lt. 10) then
         nHypoX = n2
         phypoX = real(1.0/nHypoX)
         nHypoXstep = 1
      else
         nHypoX = int(n2/10)*10
         phypoX = real(1.0/10.0)
         nHypoXstep = int(n2/10)
      endif

c     Compute the step sizes down dip
      if (n1 .lt. 10) then
         nHypoZ = n1
         phypoZ = real(1.0/nHypoZ)
         nHypoZstep = 1
      else
         nHypoZ = int(n1/10)*10
         phypoZ = real(1.0/10.0)
         nHypoZstep = int(n1/10)
      endif

      return
      end

c -------------------------------------------------------------------
      subroutine DetermDist (hAScell, hDDcell, icellRupStrike, icellRupdip,
     1                      fltgrid_x, fltgrid_y, fltgrid_z, n2, n1, dipavgd,
     2                      iLocY, iLocX, rupLen, rupWidth, x0, y0, z0,
     3                      edist, hdist, slit, azp1p2, step, dlit, phiang, FltGrid_rRup,
     4                      s2site, Rx, astrike)

      implicit none
      include 'pfrisk.h'

      integer hAScell, hDDcell, icellRupStrike, icellRupDip
      integer n2, n1, iLocx, iLocY, icell
      real  fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_x(MAXFLT_DD,MAXFLT_AS),
     1      fltGrid_y(MAXFLT_DD,MAXFLT_AS), dx, dy, dz, edist, hdist
      real dist, step, phiang, dlit, strikeX, strikeY, astrike
      real fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), Rx, dipavgd, x0, y0, z0, rupLen,
     1     rupWidth, azp1p2, slit, s2site, c1temp

C     Variable: hDDcell, hAScell --> Location of hypocenter locations down-dip and along strike.
C               iLocY, iLocX --> number of rupture areas to complete fill fault plane.
C               n2, n1 --> number of cells down-dip and along strike.
C               iCellRupdip, iCellRupStirke --> cell location down-dip and along strike for closest point to station.

C      Compute average strike for given fault plane.
       strikeX = fltgrid_x(1,n2) - fltgrid_x(1,1)
       strikeY = fltgrid_y(1,n2) - fltgrid_y(1,1)
       if (strikeX .eq. 0.0) then
          astrike = 0.0
       else
          astrike = atan2(strikeX,strikeY)
       endif

C     Compute distance along strike between closest cell on fault plane and hypocenter location.
      dx = (fltgrid_x(hDDcell,hAScell) - fltgrid_x(hDDcell,icellRupStrike))
      dy = (fltgrid_y(hDDcell,hAScell) - fltgrid_y(hDDcell,icellRupStrike))
      slit = sqrt( dx**2.0 + dy**2.0)

C     Compute the downdip distance dlit.
      dlit = real(step*(abs(icellRupDip-hDDcell)))

C     Compute the angle Phi between Hypocenter and Station location.
      phiang = (180.0/3.14159)*atan2(Rx,dlit) - (90.0 - dipavgd)

C     Compute Hypocentral and Epicentral Distances.
      dx = fltGrid_x(hDDcell,hAScell) - x0
      dy = fltGrid_y(hDDcell,hAScell) - y0
      dz = fltGrid_z(hDDcell,hAScell) - z0
      hdist = sqrt ( dx*dx + dy*dy + dz*dz)
      edist = sqrt ( dx*dx + dy*dy )

C     Compute azimuth between epicenter and station location.
      dx = x0 - fltGrid_x(hDDcell,hAScell)
      dy = y0 - fltGrid_y(hDDcell,hAScell)
      azp1p2 =  atan2(dx,dy)*(180.0/3.14159) - astrike*180/3.14159

C     Compute azimuth between closest point and site.
      dx = x0 - fltGrid_x(iCellRupdip,iCellRupstrike)
      dy = y0 - fltGrid_y(iCellRupdip,iCellRupstrike)
      s2site =  atan2(dx,dy)*(180.0/3.14159) - astrike*180/3.14159


C     Check to see if the slit value is greater than the limited c1 value. If so recompute values.
C     Case where hypocenter is further down strike than closest point cell.
      c1temp = 0.0
      if (hAScell .gt. icellRupstrike) then
         do icell=icellrupstrike, hAScell-1, 1
            dx = (fltgrid_x(icellRupdip,icell) - fltgrid_x(icellRupdip,icell+1))
            dy = (fltgrid_y(icellRupdip,icell) - fltgrid_y(icellRupdip,icell+1))
            dist = sqrt( dx**2.0 + dy**2.0)
            c1temp = c1temp + dist
         enddo
C     Case where closest point cell is further down strike than hypocenter.
      else
         do icell=hAScell,icellrupstrike-1, 1
            dx = (fltgrid_x(icellRupdip,icell) - fltgrid_x(icellRupdip,icell+1))
            dy = (fltgrid_y(icellRupdip,icell) - fltgrid_y(icellRupdip,icell+1))
            dist = sqrt( dx**2.0 + dy**2.0)
            c1temp = c1temp + dist
         enddo

      endif

      return
      end

c -------------------------------------------------------------------

      subroutine Get_plane_dist (x, y, z, x0, y0, z0, insideFlag, dist)

      implicit none

c     declarations passed in
      real x(5), y(5), z(5), x0, y0, z0

c     declarations passed out
      integer insideFlag
      real dist

c     declarations only used in this subroutine
      integer i, i1, nSeg
      real dist1

c     Determine if the site is inside the surface projection of the cell
c     boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, x, y, x0, y0, insideFlag )

c     Compute the shortest dist to each edge
      dist = 1.0e30
      do i=1,4
         i1 = i+1
         if ( i1 .gt. 4 ) i1=1
         call Calc_LineSeg_dist ( x(i), y(i), z(i), x(i1),
     1       y(i1), z(i1), x0, y0, z0, dist1 )
         if ( dist1 .lt. dist ) then
            dist = dist1
         endif
      enddo
      return
      end

c ---------------------------------------------------------------

      subroutine Calc_LineSeg_dist ( x1, y1, z1, x2, y2, z2, x0, y0,
     1           z0, dist )

      implicit none

      real x0, x1, x2, y0, y1, y2, z0, z1, z2, dist
      real t1, t2, x, y, z, L, L1, L2, d1, d2

c     Find shortest distance to line (without ends)
c     Interesection at (x,y,z)
      if ( z1 .ne. z2 ) then
         t1 = (x2-x1)/(z2-z1)
         t2 = (y2-y1)/(z2-z1)
         z =  (z0 - (-z1*t1 + x1 - x0)*t1  - (-z1*t2 + y1 - y0)*t2 )
     1     / ( t1**2 + t2**2 + 1 )
         x = t1 * (z-z1) + x1
         y = t2 * (z-z1) + y1
      elseif ( y1 .ne. y2 ) then
         z = z1
         t1 = (x2-x1)/(y2-y1)
         y = (y0 - (-y1*t1 + x1 - x0)*t1) / (t1**2 + 1)
         x = t1 * (y-y1) + x1
      else
         z = z1
         y = y1
         x = x0
      endif
      dist = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )

c     Check if intersection is outside of edge
      L = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
      L1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 )
      L2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2 )
      if ( L1 .le. L .and. L2 .le. L ) then
         return
      endif

c     Intersection does not hit segment
      d1 = sqrt( (x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2 )
      d2 = sqrt( (x0-x2)**2 + (y0-y2)**2 + (z0-z2)**2 )
      dist = min ( d1, d2 )

      return
      end

c -----------------------------

      subroutine Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0,
     1           insideFlag )

      implicit none

c     declarations passed in
      integer nSeg
      real xSeg(1), ySeg(1), x0, y0

c     declarations passed out
      integer insideFlag

c     declarations only used within subroutine
      integer i
      real twoPi, pi, theta1, theta2, dTheta, dx1, dx2, dy1, dy2,
     1     sumTheta, test, tol

c     This subroutine determines if a point (x0,y0) is inside of
c     the polygon given by xSeg, ySeg.

      pi = 3.1415926
      twoPi = 2. * pi
      sumTheta = 0.

      do i=1,nSeg

c       Compute Azimuth to ends of segments
        dy1 = ySeg(i) - y0
        dy2 = ySeg(i+1) - y0
        dx1 = xSeg(i) - x0
        dx2 = xSeg(i+1) - x0
        theta1 = atan2 ( dy1, dx1 )
        theta2 = atan2 ( dy2, dx2 )
        dTheta = theta2 - theta1

c       Check if theta range is greater than pi (wrap around)
        if ( dTheta .gt. pi ) then
           dTheta = dTheta - twoPi
        elseif ( dTheta .lt. -pi ) then
           dTheta = dTheta + twoPi
        endif

c       Compute sum of azimuth ranges
        sumTheta = sumTheta + dTheta
      enddo

c     Determine if point is inside the polygon
c     If sumTheta = +-2 pi , then point is inside
      test = abs ( abs(sumTheta) - twoPi )
      tol = 0.01

      if ( test .lt. tol ) then
         insideFlag = 1
      else
         insideFlag = 0
      endif
      return
      end

c -------------------------------------------------------------

      subroutine CalcPlaneDist ( x0, y0, z0, x, y, z, dist )

      implicit none

      integer MAXPTS, MAXTERM, arow, acol, i, nterm, npts
      parameter (MAXPTS=3, MAXTERM=3)
      real*8 A(MAXPTS,MAXTERM), b(MAXPTS,1)
      real*8 xhat(MAXPTS,1), eps, deter, work(1000)
      real x(5), y(5), z(5), x0, y0, z0, dist, cx, cy, cz, c

      eps = 1.0e-11
      arow = MAXPTS
      acol = MAXTERM

c     Compute the equation for the plane
      if ( y(1) .eq. y(2) .and. z(1) .eq. z(2) ) then
         if ( z(1) .eq. z(2) .and. z(2) .eq. z(3) ) then
             cx = 0.
             cy = 0.
             cz = 1.
             c = -z(1)
         else
             do i=1,3
                A(i,1) = x(i)
                A(i,2) = z(i)
                A(i,3) = 1.
                b(i,1) = -y(i)
            enddo
            nTerm = 3
            nPts = 3
            call simul ( nterm, A, work, eps, -1, acol, deter )
            call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1         acol )
            cy = 1.
            cx = xhat(1,1)
            cz = xhat(2,1)
            c = xhat(3,1)
         endif
      else
         do i=1,3
            A(i,1) = y(i)
            A(i,2) = z(i)
            A(i,3) = 1.
            b(i,1) = -x(i)
         enddo
      nTerm = 3
      nPts = 3
      call simul ( nterm, A, work, eps, -1, acol, deter )
      call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1       acol )
      cx = 1.
      cy = xhat(1,1)
      cz = xhat(2,1)
      c = xhat(3,1)
      endif

c     Compute distance from point to plane
      dist = abs( cx*x0 + cy*y0 + cz*z0 + c) /
     1       sqrt( cx**2 + cy**2 + cz**2 )

      return
      end

c --------------------------------------------------------------------


      subroutine Get_plane_dist2 ( xRup, yRup, zRup, iSeg,
     1           x0, y0, z0, dip0, dist, xclp, yclp, zclp )

      implicit none

      integer i, i1, nSeg, insideFlag, iSeg
      real xRup(4,1), yRup(4,1), zRup(4,1), x(4), y(4), z(4), xSeg(5),
     1     ySeg(5), dist, dip0, x0, y0, z0, dist1, pi, dip, theta,
     2     theta1, sin1, cos1, tan1, d1, d2, xclp, yclp, zclp, xclp2,
     3     yclp2, zclp2

      pi = 3.1415926

      do i=1,4
         x(i) = xRup(i,iSeg)
         y(i) = yRup(i,iSeg)
         z(i) = zRup(i,iSeg)
      enddo

c     Set angles
      dip = dip0/180 * pi
      theta = atan2 ( y(2)-y(1), x(2)-x(1) )
      theta1 = theta + pi/2.
      if ( dip .lt. 0 ) then
         theta1 = theta1 + pi
      endif
      cos1 = cos(theta1)
      sin1 = sin(theta1)
      tan1 = tan(dip)

c     Set boundary for points for which the closest point is to
C     the plane

C Note: Reorders points in clockwise rotation.
      d1 = z(1) * abs(tan1)
      d2 = z(3) * abs(tan1)
      xSeg(1) = x(1) + d1*cos1
      xSeg(2) = x(2) + d1*cos1
      xSeg(3) = x(4) + d2*cos1
      xSeg(4) = x(3) + d2*cos1
      xSeg(5) = xSeg(1)
      ySeg(1) = y(1) + d1*sin1
      ySeg(2) = y(2) + d1*sin1
      ySeg(3) = y(4) + d2*sin1
      ySeg(4) = y(3) + d2*sin1
      ySeg(5) = ySeg(1)

c     Determine if the site is inside this boundary
      nSeg = 4
      call Inside_OutSide ( nSeg, xSeg, ySeg, x0, y0, insideFlag )

c     Compute closest distance for inside
      if ( insideFlag .eq. 1 ) then
         call CalcPlaneDist2 ( x0, y0, z0, x, y, z, dist, xclp, yclp, zclp )
         return
      endif

c     Compute the shortest dist to each edge
      dist = 1.0e30
      xclp = 1.0e30
      yclp = 1.0e30
      zclp = 1.0e30
      do i=1,4
         i1 = i+1
         if ( i1 .gt. 4 ) i1=1
         call Calc_LineSeg_dist2 ( x(i), y(i), z(i), x(i1),
     1       y(i1), z(i1), x0, y0, z0, dist1, xclp2, yclp2, zclp2 )

         if ( dist1 .lt. dist ) then
            dist = dist1
            xclp = xclp2
            yclp = yclp2
            zclp = zclp2
         endif
      enddo
      return
      end

c ---------------------------------------------------------------

      subroutine Calc_LineSeg_dist2 ( x1, y1, z1, x2, y2, z2, x0, y0,
     1           z0, dist, x, y, z )

      implicit none

      real x0, x1, x2, y0, y1, y2, z0, z1, z2, dist, t1, t2, x, y, z,
     1     L, L1, L2, d1, d2

c     Find shortest distance to line (without ends)
c     Interesection at (x,y,z)
      if ( z1 .ne. z2 ) then
         t1 = (x2-x1)/(z2-z1)
         t2 = (y2-y1)/(z2-z1)
         z =  (z0 - (-z1*t1 + x1 - x0)*t1  - (-z1*t2 + y1 - y0)*t2 )
     1     / ( t1**2 + t2**2 + 1 )
         x = t1 * (z-z1) + x1
         y = t2 * (z-z1) + y1
      elseif ( y1 .ne. y2 ) then
         z = z1
         t1 = (x2-x1)/(y2-y1)
         y = (y0 - (-y1*t1 + x1 - x0)*t1) / (t1**2 + 1)
         x = t1 * (y-y1) + x1
      else
         z = z1
         y = y1
         x = x0
      endif
      dist = sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 )

c     Check if intersection is outside of edge
      L = sqrt( (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2 )
      L1 = sqrt( (x-x1)**2 + (y-y1)**2 + (z-z1)**2 )
      L2 = sqrt( (x-x2)**2 + (y-y2)**2 + (z-z2)**2 )
      if ( L1 .le. L .and. L2 .le. L ) then
         return
      endif

c     Intersection does not hit segment
      d1 = sqrt( (x0-x1)**2 + (y0-y1)**2 + (z0-z1)**2 )
      d2 = sqrt( (x0-x2)**2 + (y0-y2)**2 + (z0-z2)**2 )
      dist = min ( d1, d2 )

      return
      end

c -------------------------------------------------------------

      subroutine CalcPlaneDist2 ( x0, y0, z0, x, y, z, dist, cx, cy, cz )

      implicit none

      integer MAXPTS, MAXTERM, arow, acol, i, nterm, npts
      parameter (MAXPTS=3, MAXTERM=3)
      real x(4), y(4), z(4), x0, y0, z0, dist, cx, cy, cz, c
      real*8 A(MAXPTS,MAXTERM), b(MAXPTS,1), xhat(MAXPTS,1), eps,
     1       deter, work(1000)

      eps = 1.0e-11
      arow = MAXPTS
      acol = MAXTERM

c     Compute the equation for the plane
      if ( y(1) .eq. y(2) .and. z(1) .eq. z(2) ) then
         if ( z(1) .eq. z(2) .and. z(2) .eq. z(3) ) then
             cx = 0.
             cy = 0.
             cz = 1.
             c = -z(1)
         else
             do i=1,3
                A(i,1) = x(i)
                A(i,2) = z(i)
                A(i,3) = 1.
                b(i,1) = -y(i)
            enddo
            nTerm = 3
            nPts = 3
            call simul ( nterm, A, work, eps, -1, acol, deter )
            call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1         acol )
            cy = 1.
            cx = xhat(1,1)
            cz = xhat(2,1)
            c = xhat(3,1)
         endif
      else
         do i=1,3
            A(i,1) = y(i)
            A(i,2) = z(i)
            A(i,3) = 1.
            b(i,1) = -x(i)
         enddo
      nTerm = 3
      nPts = 3
      call simul ( nterm, A, work, eps, -1, acol, deter )
      call mult ( A, acol, nterm, npts, b, arow, npts, 1, xhat,
     1       acol )
      cx = 1.
      cy = xhat(1,1)
      cz = xhat(2,1)
      c = xhat(3,1)
      endif

c     Compute distance from point to plane
      dist = abs( cx*x0 + cy*y0 + cz*z0 + c) /
     1       sqrt( cx**2 + cy**2 + cz**2 )

      return
      end




       subroutine GC2 (iLocX, iLocY, n2, n1, fltgrid_x1, fltgrid_y1, fltgrid_z1,
     1                 fltgrid_x2, fltgrid_y2, fltgrid_z2, fltgrid_x3, fltgrid_y3,
     2                 fltgrid_z3, fltgrid_x4, fltgrid_y4, fltgrid_z4, Rx, Ry, Ry0,
     3                 HWFlag, dipavgd)

       implicit none
       include 'pfrisk.h'

c      declarations passed in
       integer iLocX, iLocY, n2, n1
       real fltgrid_x1(MAXFLT_DD,MAXFLT_AS), fltgrid_y1(MAXFLT_DD,MAXFLT_AS),
     1      fltgrid_z1(MAXFLT_DD,MAXFLT_AS), fltgrid_x2(MAXFLT_DD,MAXFLT_AS),
     2      fltgrid_y2(MAXFLT_DD,MAXFLT_AS), fltgrid_z2(MAXFLT_DD,MAXFLT_AS),
     3      fltgrid_x3(MAXFLT_DD,MAXFLT_AS), fltgrid_y3(MAXFLT_DD,MAXFLT_AS),
     4      fltgrid_z3(MAXFLT_DD,MAXFLT_AS), fltgrid_x4(MAXFLT_DD,MAXFLT_AS),
     5      fltgrid_y4(MAXFLT_DD,MAXFLT_AS), fltgrid_z4(MAXFLT_DD,MAXFLT_AS)

c      declarations passed out
       integer HWFlag
       real Rx, Ry, Ry0, dipavgd

c      declarations only used within subroutine
       integer inorm, irup, n3, iGC2, tflag, uflag, i
       real rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), rup_z(MAXFLT_AS),
     1      rup_xb(MAXFLT_AS), rup_yb(MAXFLT_AS), rup_zb(MAXFLT_AS),
     2      Seg_length(MAXFLT_AS), Strike_slope(MAXFLT_AS),
     3      Normal_slope(MAXFLT_AS), GC2_ruplength, a, b, P90_x(MAXFLT_AS),
     4      P90_y(MAXFLT_AS), Site_x, Site_y, t_local(MAXFLT_AS),
     5      u_local(MAXFLT_AS), Seg_weight(MAXFLT_AS), Seg_weight_t(MAXFLT_AS),
     6      sum_Weight, rec_Weight, Seg_x(MAXFLT_AS), Seg_wxu(MAXFLT_AS), sum_Swt,
     7      sum_Swxu, Global_T, Global_U, dipX, dipY, mdipX(MAXFLT_AS),
     8      mdipY(MAXFLT_AS), mdip(MAXFLT_AS), dipavgr

c      save rupture grid cell locations in new arrays
       inorm = 0
       do irup=iLocX, n2
         inorm = inorm + 1
         rup_x(inorm) = fltgrid_x1(iLocY,irup)
         rup_y(inorm) = fltgrid_y1(iLocY,irup)
         rup_z(inorm) = fltgrid_z1(iLocY,irup)
         rup_xb(inorm) = fltgrid_x4(n1,irup)
         rup_yb(inorm) = fltgrid_y4(n1,irup)
         rup_zb(inorm) = fltgrid_z4(n1,irup)
         if (irup .eq. n2) then
           inorm = inorm + 1
           rup_x(inorm) = fltgrid_x2(iLocY,irup)
           rup_y(inorm) = fltgrid_y2(iLocY,irup)
           rup_z(inorm) = fltgrid_z2(iLocY,irup)
           rup_xb(inorm) = fltgrid_x3(n1,irup)
           rup_yb(inorm) = fltgrid_y3(n1,irup)
           rup_zb(inorm) = fltgrid_z3(n1,irup)
         endif
       enddo
       n3 = inorm

c       calculate the length of each segment of the rupture (a segment of the
c       rupture is a grid cell)

        Site_x = 0.0
        Site_y = 0.0

        do iGC2=1, n3-1
          Seg_length(iGC2) = sqrt(((rup_x(iGC2+1)-rup_x(iGC2))**2)+((
     1                    rup_y(iGC2+1)-rup_y(iGC2))**2))
          Strike_slope(iGC2) = (rup_y(iGC2+1)-rup_y(iGC2))/(rup_x(iGC2+1)-
     1                    rup_x(iGC2))
          Normal_slope(iGC2) = (-1)/Strike_slope(iGC2)
        enddo

c       calculate total segment length (rupture length) for Ry0 adjustment
c       later

        GC2_ruplength = 0.0
        do iGC2=1, n3-1
          GC2_ruplength = GC2_ruplength + Seg_length(iGC2)
        enddo

c       calculate x and y coordinate of point where local
c       u and t form 90 degree angle

          do iGC2=1, n3-1
            a = rup_x(iGC2+1) - rup_x(iGC2)
            b = rup_y(iGC2+1) - rup_y(iGC2)
              if (a.EQ.0) then
                P90_x(iGC2) = rup_x(iGC2)
                P90_y(iGC2) = Site_y
              else if (b.EQ.0) then
                P90_x(iGC2) = Site_x
                P90_y(iGC2) = rup_y(iGC2)
              else
                P90_x(iGC2) = (((-1)*Strike_slope(iGC2)*rup_x(iGC2)) +
     1                     rup_y(iGC2) + (Normal_slope(iGC2)*Site_x) -
     2                     Site_y)/(Normal_slope(iGC2) -
     3                     Strike_slope(iGC2))
                P90_y(iGC2) = (Strike_slope(iGC2) * (P90_x(iGC2) -
     1                     rup_x(iGC2))) + rup_y(iGC2)
              endif
          enddo

c       calculate local strike-normal coordinate t

          do iGC2=1, n3-1
              call GC2_tsign (iGC2, n1, iLocX, iLocY, rup_x, rup_y,
     1                        rup_xb, rup_yb, tflag)
              t_local(iGC2) = tflag*sqrt(((P90_x(iGC2) - Site_x)**2) +
     1                        ((P90_y(iGC2) - Site_y)**2))
          enddo

c       calculate local strike-parallel coordinate u

          do iGC2=1, n3-1
              call GC2_usign (iGC2, rup_x, rup_y, uflag)
              u_local(iGC2) = uflag*sqrt(((P90_x(iGC2) - rup_x(iGC2))**2) +
     1                         ((P90_y(iGC2) - rup_y(iGC2))**2))
          enddo

c        check for special case t=0, on extension and assign
c          alternative segment weight
c        check for special case t=0, on segment and assign dummy
c          segment weight = 0, which will not be used
c        if no special case, assign Segment weight as analytical
c          solution of 1/(r^2) evaluated from 0 to Seg_length

          do iGC2=1, n3-1
            if (t_local(iGC2).EQ.0 .and. u_local(iGC2).LT.0) then
              Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) -
     1                           (1/u_local(iGC2))
            else if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GT.Seg_length(iGC2)) then
              Seg_weight(iGC2) = (1/(u_local(iGC2) - Seg_length(iGC2))) -
     1                           (1/u_local(iGC2))
            else if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GE.0 .and.
     1        u_local(iGC2).LE.Seg_length(iGC2)) then
              Seg_weight(iGC2) = 0.0
            else
              Seg_weight(iGC2) = (1/t_local(iGC2))*((ATAN(
     1                          (Seg_length(iGC2) - u_local(iGC2))/
     2                          t_local(iGC2))) - (ATAN(((-1)*
     3                          u_local(iGC2))/t_local(iGC2))))
            endif
              Seg_weight_t(iGC2) = Seg_weight(iGC2) * t_local(iGC2)
          enddo

c       calculate the reciprocal of the sum of the segment weights,
c       rec_Weight

          sum_Weight = 0.0
          do iGC2=1, n3-1
            sum_Weight = sum_Weight + Seg_weight(iGC2)
          enddo
            rec_Weight = 1/sum_Weight

c       calculate where you are on the strike of the fault for each
c       segment and segment weight * Seg_x + u_local

        Seg_x(1)=0.0
        do iGC2=2, n3-1
          Seg_x(iGC2) = Seg_x(iGC2-1) + Seg_length(iGC2-1)
        enddo

        do iGC2=1, n3-1
          Seg_wxu(iGC2) = Seg_weight(iGC2)*(Seg_x(iGC2)+u_local(iGC2))
        enddo

c       calculate Global Coordinate T
c       check for special case t=0 on segment

          sum_Swt = 0.0
          do iGC2=1, n3-1
            sum_Swt = sum_Swt + Seg_weight_t(iGC2)
            Global_T = rec_Weight*sum_Swt
          enddo
          do iGC2=1, n3-1
            if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GE.0 .and.
     1        u_local(iGC2).LE.Seg_length(iGC2)) then
              Global_T = 0.0
            endif
          enddo

c       calculate Global Coordinate U
c       check for special case t=0 on segment

          sum_Swxu = 0.0
          do iGC2=1, n3-1
            sum_Swxu = sum_Swxu + Seg_wxu(iGC2)
            Global_U = rec_Weight*sum_Swxu
          enddo
          do iGC2=1, n3-1
            if (t_local(iGC2).EQ.0 .and. u_local(iGC2).GE.0 .and.
     1        u_local(iGC2).LE.Seg_length(iGC2)) then
              Global_U = Seg_x(iGC2)+u_local(iGC2)
            endif
          enddo

c       calculate Rx, dipavgd, and assign HWFlag from Global Coordinate T
         HWFlag = 0
         dipX = fltgrid_x4(n1,iLocX) - fltgrid_x1(iLocY,iLocX)
         dipY = fltgrid_y4(n1,iLocX) - fltgrid_y1(iLocY,iLocX)
         if (dipX .eq. 0.0 .and. dipY .eq. 0.0) then
           Rx = (-1)*(abs(Global_T))
           HWFlag = 0
           dipavgd = 90.0
         else
           Rx = Global_T
           if (Global_T .LE. 0.0) then
             HWFlag = 0
           elseif (Global_T .GT. 0.0) then
             HWFlag = 1
           endif
           dipavgr = 0.0
           do i=1, n3
             mdipX(i) = rup_xb(i) - rup_x(i)
             mdipY(i) = rup_yb(i) - rup_y(i)
             mdip(i) = (atan2((rup_zb(i)-rup_z(i)),
     1                  sqrt(mdipX(i)*mdipX(i)+mdipY(i)*mdipY(i))))/n3
             dipavgr = dipavgr + mdip(i)
           enddo
           dipavgd = dipavgr*180./3.1415926
         endif

c       calculate Ry from Global Coordinate U
         if (Global_U .LT. 0) then
           Ry = abs(Global_U) + (0.5*GC2_ruplength)
         elseif (Global_U .GE. 0 .and. Global_U .LE. GC2_ruplength) then
           Ry = abs((0.5*GC2_ruplength)-Global_U)
         elseif (Global_U .GT. GC2_ruplength) then
           Ry = Global_U - (0.5*GC2_ruplength)
         endif

c       calculate Ry0 from Global Coordinate U
         if (Global_U .LT. 0) then
           Ry0 = abs(Global_U)
         elseif (Global_U .GE. 0 .and. Global_U .LE. GC2_ruplength) then
           Ry0 = 0.0
         elseif (Global_U .GT. GC2_ruplength) then
           Ry0 = Global_U - GC2_ruplength
         endif

        return
       end

c ----------------------------------------------------------------------

       subroutine GC2_tsign (iGC2, n1, iLocX, iLocY, rup_x, rup_y, rup_xb,
     1                       rup_yb, tflag)

       implicit none
       include 'pfrisk.h'

c      declarations passed in
       integer iGC2, n1, iLocX, iLocY
       real rup_x(MAXFLT_AS), rup_y(MAXFLT_AS), rup_xb(MAXFLT_AS), rup_yb(MAXFLT_AS)

c      declarations passed out
       integer tflag

c      declarations only used within subroutine
       integer tpositive, tnegative
       real strikeX, strikeY, dipX, dipY, ddip, strike, xtemp, ytemp, xtendp, ytendp,
     1      xtendn, ytendn, xtempp, ytempp, xtestp(5), ytestp(5), xtempn,
     2      ytempn, xtestn(5), ytestn(5)

c      Compute direction of strike and direction of dip (ddip)
c      For pure strike slip fault, use right hand rule for dummy dip direction
       strikeX = rup_x(iGC2+1) - rup_x(iGC2)
       strikeY = rup_y(iGC2+1) - rup_y(iGC2)
       strike = atan2(strikeY,strikeX)

       dipX = rup_xb(iGC2) - rup_x(iGC2)
       dipY = rup_yb(iGC2) - rup_y(iGC2)
         if (dipX .eq. 0.0 .and. dipY .eq. 0.0) then
           ddip = strike - 3.141592653590/2.0
         else
         ddip = atan2(dipY,dipX)
         endif

c      Extend the first point of the segment by 1000 km in each
c      of the strike directions
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       xtemp = 1000.0*cos(strike)
       ytemp = 1000.0*sin(strike)
       xtendp = rup_x(iGC2) + xtemp
       ytendp = rup_y(iGC2) + ytemp
       xtendn = rup_x(iGC2) - xtemp
       ytendn = rup_y(iGC2) - ytemp

c      Now determine if site is located on the hanging wall (t+) or on
c      the footwall (t-) relative to the segment.
c      Reset tflag for each segment.
       tpositive = 0
       tnegative = 0

c      Extend the first point of the segment by 1000 km in the direction
c      of dip. Set up testing points for t+

       xtempp = 1000*cos(ddip)
       ytempp = 1000*sin(ddip)

       xtestp(1) = xtendp
       ytestp(1) = ytendp
       xtestp(2) = xtendn
       ytestp(2) = ytendn
       xtestp(3) = xtestp(2) + xtempp
       ytestp(3) = ytestp(2) + ytempp
       xtestp(4) = xtestp(1) + xtempp
       ytestp(4) = ytestp(1) + ytempp
       xtestp(5) = xtestp(1)
       ytestp(5) = ytestp(1)

c      Extend the first point of the segment by 1000 km in the direction
c      opposite of dip. Set up testing points for t-

       xtempn = 1000*cos(ddip - 3.141592653590)
       ytempn = 1000*sin(ddip - 3.141592653590)

       xtestn(1) = xtendp
       ytestn(1) = ytendp
       xtestn(2) = xtendn
       ytestn(2) = ytendn
       xtestn(3) = xtestn(2) + xtempn
       ytestn(3) = ytestn(2) + ytempn
       xtestn(4) = xtestn(1) + xtempn
       ytestn(4) = ytestn(1) + ytempn
       xtestn(5) = xtestn(1)
       ytestn(5) = ytestn(1)

c      Check to see if site is located in the t+ testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestp, ytestp, 0.0, 0.0, tpositive)

c      Check to see if site is located in the t- testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestn, ytestn, 0.0, 0.0, tnegative)

       if (tpositive.eq.1.) then
         tflag = 1.
       elseif (tnegative.eq.1.) then
         tflag = -1.
       endif

        return
       end

c ----------------------------------------------------------------------

       subroutine GC2_usign (iGC2, rup_x, rup_y, uflag)

       implicit none
       include 'pfrisk.h'

       integer iGC2, uflag
       real rup_x(MAXFLT_AS), rup_y(MAXFLT_AS)

       integer upositive, unegative
       real strikeX, strikeY, strike, normal, xtemp, ytemp, xtendp, ytendp,
     1      xtendn, ytendn, xtempp, ytempp, xtestp(5), ytestp(5), xtempn,
     2      ytempn, xtestn(5), ytestn(5)

c      Compute strike and strike-normal for segment
       strikeX = rup_x(iGC2+1) - rup_x(iGC2)
       strikeY = rup_y(iGC2+1) - rup_y(iGC2)
       strike = atan2(strikeY,strikeX)
       normal = strike - 3.141592653590/2.0

c      Extend the first point of the segment by 1000 km in each
c      of the strike normal directions
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       xtemp = 1000.0*cos(normal)
       ytemp = 1000.0*sin(normal)
       xtendp = rup_x(iGC2) + xtemp
       ytendp = rup_y(iGC2) + ytemp
       xtendn = rup_x(iGC2) - xtemp
       ytendn = rup_y(iGC2) - ytemp

c      Now determine if site is located in the direction of rupture
c      (u+) or in the opposite direction of rupture (u-) relative to
c      the first point of this segment. Reset uflag for each segment.
       upositive = 0
       unegative = 0

c      Extend segment end point locations by 1000 km in the direction
c      of rupture. Set up testing points for u+

       xtempp = 1000*cos(strike)
       ytempp = 1000*sin(strike)

       xtestp(1) = xtendp
       ytestp(1) = ytendp
       xtestp(2) = xtendn
       ytestp(2) = ytendn
       xtestp(3) = xtestp(2) + xtempp
       ytestp(3) = ytestp(2) + ytempp
       xtestp(4) = xtestp(1) + xtempp
       ytestp(4) = ytestp(1) + ytempp
       xtestp(5) = xtestp(1)
       ytestp(5) = ytestp(1)

c      Extend segment end point locations by 1000 km in the direction
c      opposite of rupture. Set up testing points for u-

       xtempn = 1000*cos(strike - 3.141592653590)
       ytempn = 1000*sin(strike - 3.141592653590)

       xtestn(1) = xtendp
       ytestn(1) = ytendp
       xtestn(2) = xtendn
       ytestn(2) = ytendn
       xtestn(3) = xtestn(2) + xtempn
       ytestn(3) = ytestn(2) + ytempn
       xtestn(4) = xtestn(1) + xtempn
       ytestn(4) = ytestn(1) + ytempn
       xtestn(5) = xtestn(1)
       ytestn(5) = ytestn(1)

c      Check to see if site is located in the u+ testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestp, ytestp, 0.0, 0.0, upositive)

c      Check to see if site is located in the u- testing area.
c      Site is assumed to be location: (0.0, 0.0, 0.0)
       call Inside_OutSide ( 4, xtestn, ytestn, 0.0, 0.0, unegative)

       if (upositive.eq.1.) then
         uflag = 1.
       elseif (unegative.eq.1.) then
         uflag = -1.
       endif

        return
       end

c ---------------------------------------------------------------------
      subroutine ConvertCoordinates2_f (nfp, iFlt, iCoor, grid_n, sourceType, nDD,
     1                                siteX, siteY, fLat, fLong, fZ, grid_lat,
     2                                grid_long, grid_dlat, grid_dlong, nPts1, xFlt,
     3                                yFlt, zFlt, grid_x, grid_y, grid_dx, grid_dy,
     4                                x0, y0, z0)

c     This subroutine uses the haversine formula to convert latitude/longitude to
c     x/y coordinates on the basis of a spherical earth.
c     Reference: http://www.movable-type.co.uk/scripts/latlong.html

      implicit none
      include 'pfrisk.h'

      integer nfp, iFlt, iCoor, grid_n, sourceType, nDD, nPts1
      real siteX, siteY, fLat(MAX_FLT,MAX_DD,MAX_SEG), fLong(MAX_FLT,MAX_DD,MAX_SEG),
     1     fZ(MAX_FLT,MAX_DD,MAX_SEG), grid_lat(MAX_FLT,MAX_GRID),
     2     grid_long(MAX_FLT,MAX_GRID), grid_dlat(MAX_FLT), grid_dlong(MAX_FLT),
     3     xFlt(MAX_DD,MAX_SEG), yFlt(MAX_DD,MAX_SEG), zFlt(MAX_DD,MAX_SEG), grid_x(MAX_GRID),
     4     grid_y(MAX_GRID), grid_dx, grid_dy, x0, y0, z0

      integer i, iz, iPt
      real Rearth, refLat, refLong, phi1, phi2, deltaphi, lambda1, lambda2,
     1     deltalambda, a, c, dist, bearing, xscale, yscale

c     No conversion needed for SourceType 7
      if (sourceType .eq. 7) then
        goto 100
      endif

      Rearth = 6371.

c     Load npts in fault and dip into 1-D arrays
      npts1 = nfp

      if (iCoor .eq. 1.) then
c     Convert Lat, long to km (using site coordinates as ref)
        refLat = siteY
        refLong = siteX

c       Check for grid source
        if ( sourceType .eq. 3 .or. sourceType .eq. 4 ) then
          do i=1,grid_n
              phi1 = refLat*3.1415926/180.
              phi2 = grid_lat(iFlt,i)*3.1415926/180.
              deltaphi = phi2 - phi1
              lambda1 = refLong*3.1415926/180.
              lambda2 = grid_long(iFlt,i)*3.1415926/180.
              deltalambda = lambda2 - lambda1
              a = (sin(deltaphi/2.)*sin(deltaphi/2.))+(cos(phi1)*cos(phi2)
     1            *sin(deltalambda/2.)*sin(deltalambda/2.))
              c = 2*atan2(sqrt(a),sqrt(1.-a))
              dist = Rearth*c
              if (dist .eq. 0.) then
                bearing = 0.
              else
                bearing = atan2(sin(lambda2-lambda1)*cos(phi2),cos(phi1)*
     1                    sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda2-lambda1))
              endif
              grid_x(i) = dist*sin(bearing)
              grid_y(i) = dist*cos(bearing)
          enddo
c         Use flat earth scaling for grid_dx and grid_dy only
          xscale = 111.12 * cos(reflat/180*3.1415926)
          yscale = 111.12
          grid_dx = grid_dlong(iFlt)*xscale
          grid_dy = grid_dlat(iFlt)*yscale
        else

          do iz=1,nDD
            do iPt=1,nfp
              phi1 = refLat*3.1415926/180.
              phi2 = fLat(iFlt,iz,iPt)*3.1415926/180.
              deltaphi = phi2 - phi1
              lambda1 = refLong*3.1415926/180.
              lambda2 = fLong(iFlt,iz,iPt)*3.1415926/180.
              deltalambda = lambda2 - lambda1
              a = (sin(deltaphi/2.)*sin(deltaphi/2.))+(cos(phi1)*cos(phi2)
     1            *sin(deltalambda/2.)*sin(deltalambda/2.))
              c = 2*atan2(sqrt(a),sqrt(1.-a))
              dist = Rearth*c
              if (dist .eq. 0.) then
                bearing = 0.
              else
                bearing = atan2(sin(lambda2-lambda1)*cos(phi2),cos(phi1)*
     1                    sin(phi2)-sin(phi1)*cos(phi2)*cos(lambda2-lambda1))
              endif
              xFlt(iz,iPt) = dist*sin(bearing)
              yFlt(iz,iPt) = dist*cos(bearing)
              zFlt(iz,iPt) = fZ(iFlt,iz,iPt)

            enddo
          enddo
        endif

      else
c     Units are in km already
c     Adjust sources so that site is at 0,0
      refLat = siteY
      refLong = siteX

        if ( sourceType .eq. 3 .or. sourceType .eq. 4 ) then
          do i=1,grid_n
            grid_y(i) = grid_lat(iFlt,i) - reflat
            grid_x(i) = grid_long(iFlt,i) - reflong
          enddo
          grid_dx = grid_dlong(iFlt)
          grid_dy = grid_dlat(iFlt)
        else

          do iz=1,nDD
            do iPt=1,nfp
              xFlt(iz,iPt) = fLong(iFlt,iz,iPt) - reflong
              yFlt(iz,iPt) = fLat(iFlt,iz,iPt) - reflat
              zFlt(iz,iPt) = fZ(iFlt,iz,iPt)
            enddo
          enddo
        endif

      endif

      x0 = 0.
      y0 = 0.
      z0 = 0.

 100  continue

      return
      end

c ---------------------------------------------------------------------

      subroutine calcFltGrid_f ( xFlt, yFlt, zFlt, npts, nDD, fltGrid_x, fltGrid_y,
     1               fltGrid_z, nfltGrid, fltGrid_a, fltGrid_w, x0, y0, z0,
     2               fltGrid_Rrup, fltGrid_Rjb, faultArea, faultLen, faultW,
     3               step, fltGrid_fLen, fltGrid_x1, fltGrid_y1,
     4               fltGrid_z1, fltGrid_x2, fltGrid_y2, fltGrid_z2, fltGrid_x3,
     5               fltGrid_y3, fltGrid_z3, fltGrid_x4, fltGrid_y4, fltGrid_z4 )

      implicit none
      include 'pfrisk.h'

      integer nDD, n2(MAXFLT_DD), nfltGrid(2), n1(MAXFLT_AS), nPts, j, j1,
     1        j2, k, k1, i, i1, i2, ii, i0, j0, insideFlag, nx1, ny1, nn
      real fltGrid_x(MAXFLT_DD,MAXFLT_AS), fltGrid_y(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_z(MAXFLT_DD,MAXFLT_AS), fltGrid_fLen(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_w(MAXFLT_DD,MAXFLT_AS), fltGrid_a(MAXFLT_DD,MAXFLT_AS),
     3     fltGrid_Rrup(MAXFLT_DD,MAXFLT_AS), fltGrid_Rjb(MAXFLT_DD,MAXFLT_AS),
     4     faultArea, xFlt(MAX_DD,MAX_SEG), yFlt(MAX_DD,MAX_SEG)
      real zFlt(MAX_DD,MAX_SEG), x1n1, x2n1, y1n1, y2n1, z1n1, z2n1, x1n2,
     1     x2n2, y1n2, y2n2, z1n2, z2n2, dxn1, dyn1, dzn1, dxn2, dyn2,
     2     dzn2, x(5), y(5), z(5), x0, y0, z0, faultLen, dx, dy, dz, width1,
     3     dist, dist1, dist2, x1, y1, z1, x2, y2, z2, xLtop, yLtop, zLtop,
     4     xLbot, yLbot, zLbot, w, step, faultW, sum1, len1
      real fltGrid_x1(MAXFLT_DD,MAXFLT_AS), fltGrid_y1(MAXFLT_DD,MAXFLT_AS),
     1     fltGrid_z1(MAXFLT_DD,MAXFLT_AS), fltGrid_x2(MAXFLT_DD,MAXFLT_AS),
     2     fltGrid_y2(MAXFLT_DD,MAXFLT_AS), fltGrid_z2(MAXFLT_DD,MAXFLT_AS),
     3     fltGrid_x3(MAXFLT_DD,MAXFLT_AS), fltGrid_y3(MAXFLT_DD,MAXFLT_AS),
     4     fltGrid_z3(MAXFLT_DD,MAXFLT_AS), fltGrid_x4(MAXFLT_DD,MAXFLT_AS),
     5     fltGrid_y4(MAXFLT_DD,MAXFLT_AS), fltGrid_z4(MAXFLT_DD,MAXFLT_AS)

      real*8 dfaultArea

c     Find the widest part of the fault and set the number of grid points down dip
      sum1 = 0.
      nn = 0

      do i=1,nDD-1
        width1 = 0.
        do j=1,npts
          dx = xFlt(i+1,j) - xFlt(i,j)
          dy = yFlt(i+1,j) - yFlt(i,j)
          dz = zFlt(i+1,j) - zFlt(i,j)
          w = sqrt( dx**2 + dy**2 + dz**2 )
          if ( w .gt. width1 ) width1 = w
          sum1 = sum1 + real(w/npts)
          nn = nn + 1
        enddo
        n2(i) = nint( width1 / step )
        if (n2(i).eq.0) then
          n2(i) = 1
        endif
      enddo

      faultW = sum1
      nx1 = 0
      ny1 = 0
      ii = 0
      i0 = 0
      j0 = 0
      dfaultArea = 0
      faultLen = 0
      do i=1,nDD-1
        ny1 = ny1 + n2(i)

        do j=2,npts
          if ( i .eq. 1 ) then
            xLtop =  xFlt(i,j) - xFlt(i,j-1)
            yLtop =  yFlt(i,j) - yFlt(i,j-1)
            len1 = sqrt(xLtop**2 + yLtop**2)
            n1(j) = nint(len1/step)
            if (n1(j) .eq. 0) then
              n1(j) = 1
            endif
            nx1 = nx1 + n1(j)
            faultLen = faultLen + len1
          endif

          xLtop =  xFlt(i,j) - xFlt(i,j-1)
          yLtop =  yFlt(i,j) - yFlt(i,j-1)
          zLtop =  zFlt(i,j) - zFlt(i,j-1)
          xLbot =  xFlt(i+1,j) - xFlt(i+1,j-1)
          yLbot =  yFlt(i+1,j) - yFlt(i+1,j-1)
          zLbot =  zFlt(i+1,j) - zFlt(i+1,j-1)

          do j1=1,n1(j)
            j2 = j1 + j0

c           these x1-z2s are for the center point
            x1 = xFlt(i,j-1) + (XLtop/n1(j))*(j1-1+0.5)
            y1 = yFlt(i,j-1) + (yLtop/n1(j))*(j1-1+0.5)
            z1 = zFlt(i,j-1) + (zLtop/n1(j))*(j1-1+0.5)
            x2 = xFlt(i+1,j-1) + (XLbot/n1(j))*(j1-1+0.5)
            y2 = yFlt(i+1,j-1) + (yLbot/n1(j))*(j1-1+0.5)
            z2 = zFlt(i+1,j-1) + (zLbot/n1(j))*(j1-1+0.5)

            dx = (x2 - x1) / n2(i)
            dy = (y2 - y1) / n2(i)
            dz = (z2 - z1) / n2(i)

c           these x1-z2s are for node points 1 and 4
            x1n1 = xFlt(i,j-1) + (XLtop/n1(j))*(j1-1)
            x2n1 = xFlt(i+1,j-1) + (XLbot/n1(j))*(j1-1)
            y1n1 = yFlt(i,j-1) + (yLtop/n1(j))*(j1-1)
            y2n1 = yFlt(i+1,j-1) + (yLbot/n1(j))*(j1-1)
            z1n1 = zFlt(i,j-1) + (zLtop/n1(j))*(j1-1)
            z2n1 = zFlt(i+1,j-1) + (zLbot/n1(j))*(j1-1)

            dxn1 = (x2n1 - x1n1) / n2(i)
            dyn1 = (y2n1 - y1n1) / n2(i)
            dzn1 = (z2n1 - z1n1) / n2(i)

c           these x1-z2s are for node points 2 and 3
            x1n2 = xFlt(i,j-1) + (XLtop/n1(j))*(j1)
            x2n2 = xFlt(i+1,j-1) + (XLbot/n1(j))*(j1)
            y1n2 = yFlt(i,j-1) + (yLtop/n1(j))*(j1)
            y2n2 = yFlt(i+1,j-1) + (yLbot/n1(j))*(j1)
            z1n2 = zFlt(i,j-1) + (zLtop/n1(j))*(j1)
            z2n2 = zFlt(i+1,j-1) + (zLbot/n1(j))*(j1)

            dxn2 = (x2n2 - x1n2) / n2(i)
            dyn2 = (y2n2 - y1n2) / n2(i)
            dzn2 = (z2n2 - z1n2) / n2(i)

c           Set coordinates of center of grid and width and area of grid
            do i1=1,n2(i)
              i2 = i1 + i0

              if ( j2 .gt. MAXFLT_AS ) then
                write (*,'( 2x,''Error: increase dimension of MAXFLT_AS'')')
                stop 99
              endif
              if ( i2 .gt. MAXFLT_DD ) then
                write (*,'( 2x,''Error: increase dimension of MAXFLT_DD'')')
                stop 99
              endif

c             these cell coordinates are for the center point
              fltGrid_x(i2,j2) = x1 + dx*(i1-1+0.5)
              fltGrid_y(i2,j2) = y1 + dy*(i1-1+0.5)
              fltGrid_z(i2,j2) = z1 + dz*(i1-1+0.5)

c             these cell coordinates are for node point 1
              fltGrid_x1(i2,j2) = x1n1 + dxn1*(i1-1)
              fltGrid_y1(i2,j2) = y1n1 + dyn1*(i1-1)
              fltGrid_z1(i2,j2) = z1n1 + dzn1*(i1-1)

c             these cell coordinates are for node point 2
              fltGrid_x2(i2,j2) = x1n2 + dxn2*(i1-1)
              fltGrid_y2(i2,j2) = y1n2 + dyn2*(i1-1)
              fltGrid_z2(i2,j2) = z1n2 + dzn2*(i1-1)

c             these cell coordinates are for node point 3
              fltGrid_x3(i2,j2) = x1n2 + dxn2*(i1)
              fltGrid_y3(i2,j2) = y1n2 + dyn2*(i1)
              fltGrid_z3(i2,j2) = z1n2 + dzn2*(i1)

c             these cell coordinates are for node point 4
              fltGrid_x4(i2,j2) = x1n1 + dxn1*(i1)
              fltGrid_y4(i2,j2) = y1n1 + dyn1*(i1)
              fltGrid_z4(i2,j2) = z1n1 + dzn1*(i1)

c   *** this needs to be fixed.  It now assumes right angles ***
              fltGrid_w(i2,j2) = sqrt(dx**2+dy**2+dz**2)
              fltGrid_a(i2,j2) = fltGrid_w(i2,j2) * sqrt( (xLtop/n1(j))**2 + (ylTop/n1(j))**2 )
              dfaultArea = dfaultArea + dble(fltGrid_a(i2,j2))

            enddo
            i0 = ii
          enddo
          j0 = j0 + n1(j)
        enddo
        ii = ii + n2(i)
        i0 = ii
        j0 = 0
      enddo
      nfltGrid(1) = ny1
      nfltGrid(2) = nx1

c     Find the distance metrics for each cell.
      do i=1,ny1
        do j=1,nx1
          x(1) = fltgrid_x1(i,j)
          y(1) = fltgrid_y1(i,j)
          z(1) = fltgrid_z1(i,j)
          x(2) = fltgrid_x2(i,j)
          y(2) = fltgrid_y2(i,j)
          z(2) = fltgrid_z2(i,j)
          x(3) = fltgrid_x3(i,j)
          y(3) = fltgrid_y3(i,j)
          z(3) = fltgrid_z3(i,j)
          x(4) = fltgrid_x4(i,j)
          y(4) = fltgrid_y4(i,j)
          z(4) = fltgrid_z4(i,j)
          x(5) = x(1)
          y(5) = y(1)
          z(5) = z(1)

c           Compute Rrup
            call Get_plane_dist (x, y, z, x0, y0, z0, insideFlag, dist)
            fltgrid_Rrup(i,j) = dist

c           Compute Rjb
            if ( insideFlag .eq. 1 ) then
              fltGrid_Rjb(i,j) = 0.0
            else
              dist2 = 1.0e30
              do k=1,4
                k1 = k+1
                if ( k1 .gt. 4 ) k1=1
                call Calc_LineSeg_dist ( x(k), y(k), 0.0, x(k1),
     1          y(k1), 0.0, x0, y0, z0, dist1 )
                if ( dist1 .lt. dist2 ) then
                  dist2 = dist1
                endif
              enddo
              fltGrid_Rjb(i,j) = dist2
            endif

        enddo
      enddo

C     Compute the fault lengths for each cell relative to start of fault.
      do i=1,ny1
        dx = fltGrid_x(i,1) - fltGrid_x(i,2)
        dy = fltGrid_y(i,1) - fltGrid_y(i,2)
        fltGrid_fLen(i,1) = sqrt( dx**2 + dy**2 )
        do j=2,nx1
          dx = fltGrid_x(i,j-1) - fltGrid_x(i,j)
          dy = fltGrid_y(i,j-1) - fltGrid_y(i,j)
          fltGrid_fLen(i,j) = fltGrid_fLen(i,j-1) + sqrt( dx**2 + dy**2 )
        enddo
      enddo

      faultArea = real(dfaultArea)

      return
      end

c -------------------------------------------------------------------

      subroutine SetFltBottom_f (iCoor, iFlt, nfp, dip2, faultThick, fZ,
     1                         flat, flong, nDD)

      implicit none
      include 'pfrisk.h'

c      declarations passed in
       integer iCoor, iFlt, nfp(MAX_FLT)
       real dip2, faultThick

c      declarations passed out
       integer nDD(MAX_FLT)

c      declarations passed in and out (changed)
       real fZ(MAX_FLT,MAX_DD,MAX_SEG), fLat(MAX_FLT,MAX_DD,MAX_SEG),
     1      fLong(MAX_FLT,MAX_DD,MAX_SEG)

c      declarations only used within subroutine
       integer ipt
       real sin_theta, top, bottom, lat1, dx, dy, strike, az1, R1, x1, y1

C     Set bottom for case in which fault is defined in latitude and longitude.
      if (iCoor .eq. 1) then
         sin_theta = sin( abs(dip2) *3.1415926/180.)
         nDD(iFlt) = 2.
         top = fZ(iFlt,1,1)
         bottom = top + faultThick
         lat1 = ( flat(iFlt,1,nfp(iflt)) + flat(iFlt,1,1) )/2. * 3.14159/180.
         dx = ( flong(iFlt,1,nfp(iflt)) - flong(iFlt,1,1) ) * cos(lat1) * 111.12
         dy = (flat(iFlt,1,nfp(iflt)) - flat(iFlt,1,1)) * 111.12
         strike = atan2(dx,dy)
         if ( dip2 .gt. 0 ) then
           az1 = strike + 3.1415926/2.
         else
           az1 = strike - 3.1415926/2.
         endif
         R1 = faultThick / tan(abs(dip2)*3.1415926/180.)
         do ipt=1,nfp(iflt)
           x1 = R1 * sin(az1)
           y1 = R1 * cos(az1)
           flong(iFlt,2,ipt) = flong(iFlt,1,ipt) + x1/(111.12*cos(lat1))
           flat(iFlt,2,ipt)  = flat(iFlt,1,ipt) + y1/111.12
           fZ(iFlt,2,ipt) = bottom
         enddo
C     Case where fault is defined in terms of X,Y coordinates in km.
      elseif (iCoor .eq. 0) then
         sin_theta = sin( abs(dip2) *3.1415926/180.)
         nDD(iFlt) = 2.
         top = fZ(iFlt,1,1)
         bottom = top + faultThick
         lat1 = 0.0
         dx = ( flong(iFlt,1,nfp(iflt)) - flong(iFlt,1,1) )
         dy = (flat(iFlt,1,nfp(iflt)) - flat(iFlt,1,1))
         strike = atan2(dx,dy)
         if ( dip2 .gt. 0 ) then
           az1 = strike + 3.1415926/2.
         else
           az1 = strike - 3.1415926/2.
         endif
         R1 = faultThick / tan(abs(dip2)*3.1415926/180.)
         do ipt=1,nfp(iflt)
           x1 = R1 * sin(az1)
           y1 = R1 * cos(az1)
           flong(iFlt,2,ipt) = flong(iFlt,1,ipt) + x1
           flat(iFlt,2,ipt)  = flat(iFlt,1,ipt) + y1
           fZ(iFlt,2,ipt) = bottom
         enddo

      endif

      return
      end

c -----------------------------------------------------------------

      subroutine mult (a,arow,m1,n1,b,brow,m2,n2,c,crow)
c     This subroutine multiplies two matrices.
c       C = A B
c     where A is m1 x n1  and  B is m2 x n2

      implicit none

      INTEGER arow, brow, crow, m1, m2, n1, n2, i, j, ii
      double precision a(arow,1), b(brow,1), c(crow,1)

      if (n1 .ne. m2) stop 99
      do 100 i=1,m1
        do 90 j=1,n2
          c(i,j) = 0.0
          do 80 ii=1,n1
            c(i,j) = c(i,j) + a(i,ii)*b(ii,j)
  80      continue
  90    continue
  100 continue
      return
      end

c-------------------------------------------------------------
      SUBROUTINE simul( N, A, X, EPS, INDIC, NRC , DETER )
C
C     WHEN INDIC IS NEGATIVE, SIMUL COMPUTES THE INVERSE OF THE N BY
C     N MATRIX A IN PLACE.  WHEN INDIC IS ZERO, SIMUL COMPUTES THE
C     N SOLUTIONS X(1)...X(N) CORRESPONDING TO THE SET OF LINEAR
C     EQUATIONS WITH AUGMENTED MATRIX OF COEFFICIENTS IN THE N BY
C     N+1 ARRAY A AND IN ADDITION COMPUTES THE INVERSE OF THE
C     COEFFICIENT MATRIX IN PLACE AS ABOVE.  IF INDIC IS POSITIVE,
C     THE SET OF LINEAR EQUATIONS IS SOLVED BUT THE INVERSE IS NOT
C     COMPUTED IN PLACE.  THE GAUSS-JORDAN COMPLETE ELIMINATION METHOD
C     IS EMPLOYED WITH THE MAXIMUM PIVOT STRATEGY.  ROW AND COLUMN
C     SUBSCRIPTS OF SUCCESSIVE PIVOT ELEMENTS ARE SAVED IN ORDER IN
C     THE IROW AND JCOL ARRAYS RESPECTIVELY.  K IS THE PIVOT COUNTER,
C     PIVOT THE ALGEBRAIC VALUE OF THE PIVOT ELEMENT, MAX
C     THE NUMBER OF COLUMNS IN A AND DETER THE DETERMINANT OF THE
C     COEFFICIENTS MATRIX.  THE SOLUTIONS ARE COMPUTED IN THE (N+1)TH
C     COLUMN OF A AND THEN UNSCRAMBLED AND PUT IN PROPER ORDER IN
C     X(1)...X(N) USING THE PIVOT SUBSCRIPT INFORMATION AVAILABLE
C     IN THE IROW AND JCOL ARRAYS.  THE SIGN OF THE DETERMINANT IS
C     ADJUSTED, IF NECESSARY, BY DETERMINING IF AN ODD OR EVEN NUMBER
C     OF PAIRWISE INTERCHANGES IS REQUIRED TO PUT THE ELEMENTS OF THE
C     JORD ARRAY IN ASCENDING SEQUENCE WHERE JORD(IROW(I)) = JCOL(I).
C     IF THE INVERSE IS REQUIRED, IT IS UNSCRAMBLED IN PLACE USING
C     Y(1)...Y(N) AS TEMPORARY STORAGE.  THE VALUE OF THE DETERMINANT
C     IS RETURNED AS THE VALUE OF THE FUNCTION.  SHOULD THE POTENTIAL
C     PIVOT OF LARGEST MAGNITUDE BE SMALLER IN MAGNITUDE THAN EPS,
C     THE MATRIX IS CONSIDERED TO BE SINGULAR AND A TRUE ZERO IS
C     RETURNED AS THE VALUE OF THE FUNCTION.
C
      IMPLICIT REAL*8(A-H, O-Z)
      real*8 A(nrc,nrc), X(N), Y(200), eps, deter
      integer IROW(200), JCOL(200), JORD(200)
C
      MAX = N
      IF ( INDIC.GE.0 )  MAX = N + 1
C
C     ..... IS N LARGER THAN 78 .....
      IF ( N.LE.78 )  GO TO 5
      WRITE (6,200)
      DETER = 0.
      RETURN
C
C     ..... BEGIN ELIMINATION PROCEDURE .....
 5    DETER = 1.
      DO 18 K = 1, N
      KM1 = K - 1
C     ..... SEARCH FOR THE PIVOT ELEMENT .....
      PIVOT = 0.
      DO 11 I = 1, N
      DO 11 J = 1, N
C     ..... SCAN IROW AND JCOL ARRAYS FOR INVALID PIVOT SUBSCRIPTS ....
      IF ( K.EQ.1 ) GO TO 9
      DO 8 ISCAN = 1, KM1
      DO 8 JSCAN = 1, KM1
      IF ( I.EQ.IROW(ISCAN) ) GO TO 11
      IF ( J.EQ.JCOL(JSCAN) ) GO TO 11
 8    CONTINUE
 9    CONTINUE
      PIVOT = A(I,J)
      IROW(K) = I
      JCOL(K) = J
 11   CONTINUE
C
C     ..... INSURE THAT SELECTED PIVOT IS LARGER THAN EPS .....
      IF ( DABS(PIVOT).GT.EPS ) GO TO 13
      DETER = 0.
      RETURN
C
C     ..... UPDATE THE DETERMINANT VALUE .....
 13   IROWK = IROW(K)
      JCOLK = JCOL(K)
      DETER = DETER*PIVOT
C
C     ..... NORMALIZE PIVOT ROW ELEMENTS .....
      DO 14 J = 1, MAX
 14   A(IROWK,J) = A(IROWK,J)/PIVOT
C
C     ..... CARRY OUT ELIMINATION AND  DEVELOP INVERSE .....
      A(IROWK,JCOLK) = 1./PIVOT
      DO 18 I = 1, N
      AIJCK = A(I,JCOLK)
      IF ( I.EQ.IROWK ) GO TO 18
      A(I,JCOLK) = - AIJCK/PIVOT
      DO 17 J = 1, MAX
 17   IF ( J.NE.JCOLK ) A(I,J) = A(I,J) - AIJCK*A(IROWK,J)
 18   CONTINUE
C
C     ..... ORDER SOLUTION VALUES (IF ANY) AND CREATE JORD ARRAY .....
      DO 20 I = 1, N
      IROWI = IROW(I)
      JCOLI = JCOL(I)
      JORD(IROWI) = JCOLI
 20   IF ( INDIC.GE.0 ) X(JCOLI) = A(IROWI,MAX)
C
C     ..... ADJUST SIGN OF DETERMINANT .....
      INTCH = 0
      NM1 = N - 1
      DO 22 I = 1, NM1
      IP1 = I + 1
      DO 22 J = IP1,N
      IF ( JORD(J).GE.JORD(I) ) GO TO 22
      JTEMP = JORD(J)
      JORD(J) = JORD(I)
      JORD(I) = JTEMP
      INTCH = INTCH + 1
 22   CONTINUE
      IF( INTCH/2*2.NE.INTCH ) DETER = - DETER
C
C     ..... IF INDIC IS POSITIVE RETURN WITH RESULTS .....
      IF ( INDIC.LE.0 ) GO TO 26
c     DETER = DETER
      RETURN
C
C     ..... IF INDIC IS NEGATIVE OR ZERO, UNSCRAMBLE THE INVERSE
C           FIRST BY ROWS .....
 26   DO 28 J = 1, N
      DO 27 I = 1, N
      IROWI = IROW(I)
      JCOLI = JCOL(I)
 27   Y(JCOLI) = A(IROWI,J)
      DO 28 I = 1, N
 28   A(I,J) = Y(I)
C     ..... THEN BY COLUMNS .....
      DO 30 I = 1, N
      DO 29 J = 1, N
      IROWJ = IROW(J)
      JCOLJ = JCOL(J)
 29   Y(IROWJ) = A(I,JCOLJ)
      DO 30 J = 1, N
 30   A(I,J) = Y(J)
C
C     ..... RETURN FOR INDIC NEGATIVE OR ZERO .....
c     DETER = DETER
      RETURN
C
C     ..... FORMAT FOR OUTPUT STATEMENT .....
 200  FORMAT( 10H0N TOO BIG )
C
      END

c----------------------------------------------------------------------
