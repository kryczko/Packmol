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
c Subroutine that computes the analytical derivatives
c
c
      subroutine geasy(n, x, g) 

      implicit none
      include 'sizes.i'
      include 'molpa.i'

      integer n
      integer idatom, iatom, irest 
      integer i, j, k, ilubar, ilugan, icart, itype, imol
      integer imoltmp
      integer iboxx, iboxy, iboxz
      integer k1, k2
      integer iratcount
      integer tot_boxes,min_box,max_box,box_counter,nl
      integer my_rank, num_proc, OMP_GET_NUM_THREADS,
     +   OMP_GET_THREAD_NUM, CHUNK
      integer t1,t2,rate
      integer mini,maxi
      double precision xc1,xc2,xc3

      double precision x(n), g(n) 
      double precision gxcar(maxatom, 3), gxctmp(maxatom, 3)
      double precision dv1beta(3), dv1gama(3), dv1teta(3), 
     +                 dv2beta(3), dv2gama(3), dv2teta(3), 
     +                 dv3beta(3), dv3gama(3), dv3teta(3) 
      double precision v1(3), v2(3), v3(3)
      double precision xbar, ybar, zbar
      double precision clength, a1, a2, a3, b1, b2, b3, 
     +                 c1, c2, d, w
      double precision xmin, ymin, zmin, xmax, ymax, zmax
      double precision xtemp, ytemp, ztemp
      double precision beta, gama, teta, cb, sb, cg, sg, ct, st
      double precision ktemp,itemp

c Reset gradients
c      call system_clock(t1, rate)

!$OMP PARALLEL PRIVATE(i) DEFAULT(SHARED)
!$OMP DO SCHEDULE (STATIC)
      do i=1,ntotat
         gxcar (i, 1) = 0.d0 
         gxcar (i, 2) = 0.d0 
         gxcar (i, 3) = 0.d0 
      end do 
!$OMP END DO
!$OMP END PARALLEL
c      call system_clock(t2)
c      gtime = gtime + real(t2-t1)/rate

c Reset boxes

      call resetboxes(nboxes,latomfirst,latomfix)

c Transform baricenter and angles into cartesian coordinates 
c Computes cartesian coordinates from vector x and coor 
 
      do itype = 1, ntype 

        if(comptype(itype)) then

c           call system_clock(t1,rate)

!$OMP PARALLEL PRIVATE(imol,xbar,ybar,zbar,beta,
!$OMP+gama,teta,v1,v2,v3,idatom,iatom,icart,my_rank,
!$OMP+iratcount,irest,clength,xmin,ymin,zmin,ilugan,ilubar,
!$OMP+xmax,ymax,zmax,a1,a2,a3,b1,b2,b3,c1,c2,d,w) DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)
        do imol = 1, nmols(itype) 

          call compiluganbar(itype,imol,ilugan,ilubar) 

          xbar = x(ilubar + 1) 
          ybar = x(ilubar + 2) 
          zbar = x(ilubar + 3) 
 
c Compute the rotation matrix 

          beta = x(ilugan + 1)
          gama = x(ilugan + 2)
          teta = x(ilugan + 3)

          call eulerrmat(beta,gama,teta,v1,v2,v3)  
    
          idatom = idfirst(itype) - 1
          do iatom = 1, natoms(itype) 
    
            call compicart(itype,imol,iatom,icart)
            idatom = idatom + 1

            call compcart(icart,xcart,
     +                    xbar,ybar,zbar,
     +                    coor(idatom,1),coor(idatom,2),coor(idatom,3),
     +                    v1,v2,v3)


c Gradient relative to the wall distace

            do iratcount = 1, nratom(icart)
               irest = iratom(icart,iratcount)
               call gwalls(icart,irest,gxcar)
            end do
          end do
        end do
!$OMP END DO

!$OMP END PARALLEL
c        call system_clock(t2)
c        gtime = gtime + real(t2-t1)/rate
        end if

c Add the atoms of this type of molecule to the lists

        if(comptype(itype)) then
           do imol = 1, nmols(itype) 

              idatom = idfirst(itype) - 1
              do iatom = 1, natoms(itype) 

                 if(.not.init1) then    
                    call compicart(itype,imol,iatom,icart)

                    xtemp = xcart(icart,1) - sizemin(1)
                    ytemp = xcart(icart,2) - sizemin(2) 
                    ztemp = xcart(icart,3) - sizemin(3) 
                    
                    iboxx = int(xtemp/boxl(1)) + 1
                    iboxy = int(ytemp/boxl(2)) + 1
                    iboxz = int(ztemp/boxl(3)) + 1
                    
                    if(xtemp.le.0) iboxx = 1
                    if(ytemp.le.0) iboxy = 1
                    if(ztemp.le.0) iboxz = 1 
                    if(iboxx.gt.nboxes(1)) iboxx = nboxes(1)
                    if(iboxy.gt.nboxes(2)) iboxy = nboxes(2)
                    if(iboxz.gt.nboxes(3)) iboxz = nboxes(3)
                    
                    latomnext(icart) = latomfirst(iboxx,iboxy,iboxz)
                    latomfirst(iboxx,iboxy,iboxz) = icart
                    
                    ibtype(icart) = itype
                    ibmol(icart) = imol  
                 end if
              end do
           end do
        end if
      end do

      if( .not. init1) then

c
c Gradient relative to minimum distance
c

      tot_boxes = nboxes(1)*nboxes(2)*nboxes(3)

c      call system_clock(t1,rate)

!$OMP PARALLEL PRIVATE(my_rank,num_proc,min_box,max_box,i,j,k,icart,
!$OMP+box_counter,ktemp,nl,itemp) DEFAULT(SHARED)

      my_rank = OMP_GET_THREAD_NUM()
      num_proc = OMP_GET_NUM_THREADS()

!$OMP DO SCHEDULE(STATIC)

      do box_counter = 1,tot_boxes
         ktemp = dble(box_counter)/dble(nboxes(1)*nboxes(2))
         call ceil(ktemp,k)
         nl = box_counter - (k-1)*nboxes(1)*nboxes(2)
         itemp = dble(nl)/dble(nboxes(2))
         call ceil(itemp,i)
         j = nl - (i-1)*nboxes(2)

         icart = latomfirst(i,j,k)
         do while ( icart .ne. 0 ) 
            
            if(comptype(ibtype(icart))) then
               
c           Interactions inside box

                call gparc(icart,latomnext(icart),gxcar)

c           Interactions of boxes that share faces

                call gparc(icart,latomfirst(i+1,j,k),gxcar)
                call gparc(icart,latomfirst(i,j+1,k),gxcar)
                call gparc(icart,latomfirst(i,j,k+1),gxcar)

c           Interactions of boxes that share axes

                call gparc(icart,latomfirst(i+1,j+1,k),gxcar)
                call gparc(icart,latomfirst(i+1,j,k+1),gxcar)
                call gparc(icart,latomfirst(i+1,j-1,k),gxcar)
                call gparc(icart,latomfirst(i+1,j,k-1),gxcar)
                call gparc(icart,latomfirst(i,j+1,k+1),gxcar)
                call gparc(icart,latomfirst(i,j+1,k-1),gxcar)

c           Interactions of boxes that share vertices

                call gparc(icart,latomfirst(i+1,j+1,k+1),gxcar)
                call gparc(icart,latomfirst(i+1,j+1,k-1),gxcar)
                call gparc(icart,latomfirst(i+1,j-1,k+1),gxcar)
                call gparc(icart,latomfirst(i+1,j-1,k-1),gxcar)

             end if

             icart = latomnext(icart)

         end do
      end do

!$OMP END DO

!$OMP END PARALLEL
c      call system_clock(t2)
c      gtime = gtime + real(t2-t1)/rate
      
      end if

c Computing the gradient using chain rule 

c      call system_clock(t1,rate)
!$OMP PARALLEL PRIVATE(i) DEFAULT(SHARED)
!$OMP DO SCHEDULE(STATIC)
      do i = 1,n
        g(i) = 0.d0 
      end do 
!$OMP END DO
!$OMP END PARALLEL
c      call system_clock(t2)
c      gtime = gtime + real(t2-t1)/rate

      do itype = 1, ntype 

        if(comptype(itype)) then

c           call system_clock(t1,rate)
   
!$OMP PARALLEL PRIVATE(imol,cb,sb,cg,sg,ct,st,beta,
!$OMP+gama,teta,idatom,iatom,icart,my_rank,dv1beta,
!$OMP+dv2beta,dv3beta,dv1gama,dv2gama,dv3gama,
!$OMP+dv1teta,dv2teta,dv3teta,k,
!$OMP+imoltmp,k1,k2) DEFAULT(SHARED)

!$OMP DO SCHEDULE(STATIC)

        do imol = 1, nmols(itype)
           
          call compiluganbar(itype,imol,k2,k1) 

          beta = x(k2 + 1) 
          gama = x(k2 + 2) 
          teta = x(k2 + 3) 

          cb = dcos(beta) 
          sb = dsin(beta) 
          cg = dcos(gama) 
          sg = dsin(gama) 
          ct = dcos(teta) 
          st = dsin(teta) 
     
          dv1beta(1) = - cb * sg * ct - sb * cg 
          dv2beta(1) = - sb * sg * ct + cb * cg 
          dv3beta(1) = 0.d0 
     
          dv1gama(1) = - sb * cg * ct - cb * sg 
          dv2gama(1) =   cb * cg * ct - sb * sg 
          dv3gama(1) =   cg * st 
     
          dv1teta(1) =   sb * sg * st 
          dv2teta(1) = - cb * sg * st 
          dv3teta(1) =   sg * ct 
           
          dv1beta(2) = - cb * cg * ct + sb * sg 
          dv2beta(2) = - sb * cg * ct - cb * sg 
          dv3beta(2) = 0.d0 
     
          dv1gama(2) =   sb * sg * ct - cb * cg 
          dv2gama(2) = - sg * cb * ct - cg * sb 
          dv3gama(2) = - sg * st 
     
          dv1teta(2) =   sb * cg * st 
          dv2teta(2) = - cb * cg * st 
    
          dv3teta(2) =   cg * ct 
     
          dv1beta(3) =   cb * st 
          dv2beta(3) =   sb * st 
          dv3beta(3) = 0.d0 
     
          dv1gama(3) = 0.d0 
          dv2gama(3) = 0.d0 
          dv3gama(3) = 0.d0 
     
          dv1teta(3) =   sb * ct 
          dv2teta(3) = - cb * ct 
          dv3teta(3) = - st 

          idatom = idfirst(itype) - 1
          do iatom = 1, natoms(itype)
          
            call compicart(itype,imol,iatom,icart)

            idatom = idatom + 1 

            do k = 1, 3 
              g(k1+k) = g(k1+k) + gxcar(icart, k) 
            end do 
     
            do k = 1, 3 
              g(k2 + 1) = g(k2 + 1) 
     +                    + (coor(idatom,1) * dv1beta(k)  
     +                    + coor(idatom, 2) * dv2beta(k)  
     +                    + coor(idatom, 3) * dv3beta(k)) 
     +                    * gxcar(icart, k) 
          
              g(k2 + 2) = g(k2 + 2) 
     +                    + (coor(idatom,1)  * dv1gama(k)  
     +                    + coor(idatom, 2) * dv2gama(k)  
     +                    + coor(idatom, 3) * dv3gama(k)) 
     +                    * gxcar(icart, k) 
          
              g(k2 + 3) = g(k2 + 3) 
     +                    + (coor(idatom,1)  * dv1teta(k)  
     +                    + coor(idatom, 2) * dv2teta(k) 
     +                    + coor(idatom, 3) * dv3teta(k)) 
     +                    * gxcar(icart, k) 
            end do           

          end do 
        end do 
!$OMP END DO

!$OMP END PARALLEL
c        call system_clock(t2)
c        gtime = gtime + real(t2-t1)/rate
        end if
      end do 

      return 

      end 

c
c Compute the main part of the gradient
c

      subroutine gparc(icart,firstjcart,gxcar2)

      implicit none

C     COMMON BLOCK
      include 'sizes.i'
      include 'molpa.i'

C     SCALAR ARGUMENTS
      integer icart,firstjcart

C     ARRAY ARGUMENTS
      double precision gxcar2(maxatom,3)

C     LOCAL SCALARS
      integer jcart
      double precision a1,a2,a3,datom,dtemp,
     +     xdiff1, xdiff2, xdiff3

      jcart = firstjcart

      do while ( jcart .ne. 0 )
        if(comptype(ibtype(jcart))) then
          if(ibmol(icart).ne.ibmol(jcart).or.
     +      ibtype(icart).ne.ibtype(jcart)) then
            a1 = xcart(icart, 1)-xcart(jcart, 1) 
            a1 = a1 * a1
            if(a1.lt.dism2) then
              a2 = xcart(icart, 2)-xcart(jcart, 2) 
              a2 = a1 + a2 * a2
              if(a2.lt.dism2) then
                a3 = xcart(icart, 3)-xcart(jcart, 3)
                datom = a2 + a3 * a3 
                if(datom.lt.dism2) then 
                  dtemp = 4.d0 * (datom - dism2)
                  xdiff1 = dtemp*(xcart(icart,1) - xcart(jcart,1)) 
                  xdiff2 = dtemp*(xcart(icart,2) - xcart(jcart,2)) 
                  xdiff3 = dtemp*(xcart(icart,3) - xcart(jcart,3)) 
!$OMP CRITICAL
                  gxcar2(icart,1)= gxcar2(icart,1) + xdiff1
                  gxcar2(jcart,1)= gxcar2(jcart,1) - xdiff1

                  gxcar2(icart,2)= gxcar2(icart,2) + xdiff2
                  gxcar2(jcart,2)= gxcar2(jcart,2) - xdiff2

                  gxcar2(icart,3)= gxcar2(icart,3) + xdiff3
                  gxcar2(jcart,3)= gxcar2(jcart,3) - xdiff3
!$OMP END CRITICAL
                end if
              end if
            end if 
          end if
        end if
        jcart = latomnext(jcart)
      end do

      return
      end
         

