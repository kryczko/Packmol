c  
c  Written by Leandro Martínez, 2009-2011.
c  Copyright (c) 2009-2011, Leandro Martínez, Jose Mario Martinez,
c  Ernesto G. Birgin.
c  
c  This program is free software; you can redistribute it and/or
c  modify it under the terms of the GNU General Public License
c  as published by the Free Software Foundation; either version 2
c  of the License, or (at your option) any later version.
c  
c
c Subroutine eulerrmat: Computes the rotation matrix from the
c                       Euler angles
c 

c Note that:
c In this routine, beta is a rotation about the y-axis
c                  gama is a rotation about the z-axis
c                  teta is a rotation about the x-axis

      subroutine eulerrmat(beta,gama,teta,v1,v2,v3)

      implicit none
      double precision beta, gama, teta
      double precision cb, sb, cg, sg, ct, st
      double precision v1(3), v2(3), v3(3)

      cb = dcos(beta) 
      sb = dsin(beta) 
      cg = dcos(gama) 
      sg = dsin(gama) 
      ct = dcos(teta) 
      st = dsin(teta)

      v1(1)=-sb * sg * ct + cb * cg 
      v1(2)=-sb * cg * ct - cb * sg 
      v1(3)= sb * st 
   
      v2(1)= cb * sg * ct + sb * cg 
      v2(2)= cb * cg * ct - sb * sg 
      v2(3)=-cb * st 

      v3(1)= sg * st 
      v3(2)= cg * st 
      v3(3)= ct   

      return
      end

c
c Subroutine compcart: Compute cartesian coordinates using
c                      the center of mass, the canonical coordinates
c                      and the rotation matrix
c      

      subroutine compcart(icart,xcart,
     +                    xbar,ybar,zbar,
     +                    xcoor,ycoor,zcoor,
     +                    v1,v2,v3)


      implicit none
      include 'sizes.i'
      integer icart
      double precision xcart(maxatom,3)
      double precision xbar, ybar, zbar
      double precision xcoor, ycoor, zcoor
      double precision v1(3), v2(3), v3(3)

      xcart(icart,1) = xbar + xcoor*v1(1) + ycoor*v2(1) + zcoor*v3(1)    
      xcart(icart,2) = ybar + xcoor*v1(2) + ycoor*v2(2) + zcoor*v3(2)    
      xcart(icart,3) = zbar + xcoor*v1(3) + ycoor*v2(3) + zcoor*v3(3)    

      return
      end

c
c Subroutine eulerfixed: This routine was added because it defines 
c                        the rotation in the "human" way, an is thus used
c                        to set the position of the fixed molecules. 
c     That means: beta is a counterclockwise rotation around x axis.
c                 gama is a counterclockwise rotation around y axis.
c                 teta is a counterclockwise rotation around z axis.
c     The other routine should better do this as well, but then we need to change
c     all the derivative calculations, just for the sake of human interpretation
c     of the rotation which, in that case, is not really important. Maybe some day.
c 

      subroutine eulerfixed(beta,gama,teta,v1,v2,v3)

      implicit none
      double precision beta, gama, teta
      double precision c1, s1, c2, s2, c3, s3
      double precision v1(3), v2(3), v3(3)

      c1 = dcos(beta) 
      s1 = dsin(beta) 
      c2 = dcos(gama) 
      s2 = dsin(gama) 
      c3 = dcos(teta) 
      s3 = dsin(teta)

      v1(1) = c2*c3
      v1(2) = c1*s3 + c3*s1*s2
      v1(3) = s1*s3 - c1*c3*s2
      v2(1) = -c2*s3
      v2(2) = c1*c3 - s1*s2*s3
      v2(3) = c1*s2*s3 + c3*s1
      v3(1) = s2
      v3(2) = -c2*s1
      v3(3) = c1*c2         

      return
      end

