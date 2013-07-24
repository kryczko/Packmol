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
c Subroutine that computes the function value
c
c
c

      subroutine feasy(x,f) 
      
      implicit none
      include 'sizes.i'
      include 'molpa.i'

      double precision v1(3), v2(3), v3(3) 
      double precision x(nn)
      double precision f,fparc,fplus,farray(nbp**3)
      double precision xtemp, ytemp, ztemp
      double precision xbar, ybar, zbar
      double precision beta, gama, teta
      double precision ktemp,itemp, ftemp
      double precision fresttmp, fdtmp, fpen(maxatom)
      double precision flast
      integer t1,t2,rate
      integer apagar(maxatom)

      integer i, j, k, l
      integer ilugan, ilubar, icart, itype, imol, iatom, idatom, 
     +        iboxx, iboxy, iboxz, ifpen
      integer tot_boxes,min_box,max_box,box_counter,nl
      integer my_rank, num_proc, OMP_GET_THREAD_NUM

      if (ntotmol .le. 4 ) then
         call feasyseq(x,f)
         return
      end if

c Reset function value

      f = 0.d0 
      frest = 0.d0
c      call system_clock(t1,rate)
!$OMP PARALLEL PRIVATE(imol) DEFAULT (SHARED)
!$OMP DO SCHEDULE(GUIDED)
      do imol=1,ntotmol
         fpen(imol) = 0.d0
      end do
!$OMP END DO
!$OMP END  PARALLEL
c      call system_clock(t2)
c      ftime = ftime + real(t2-t1)/rate

c Reset boxes

      call resetboxes(nboxes,latomfirst,latomfix)

c Transform baricenter and angles into cartesian coordinates 
c Computes cartesian coordinates from vector x and coor 

      ilubar = 0 
      ilugan = ntotmol*3 
      icart = natfix


c      call system_clock(t1,rate)

!$OMP PARALLEL PRIVATE(itype,imol,xbar,ybar,zbar,beta,gama,
!$OMP+teta,v1,v2,v3,idatom,iatom,icart,fplus,ilugan,ilubar,
!$OMP+fresttmp,ifpen) DEFAULT(SHARED)
      do itype = 1, ntype 
        if(comptype(itype)) then
           fresttmp = frest
!$OMP DO SCHEDULE(GUIDED) 
           do imol = 1, nmols(itype)
              call compiluganbar(itype,imol,ilugan,ilubar) 
              
              xbar = x(ilubar + 1) 
              ybar = x(ilubar + 2) 
              zbar = x(ilubar + 3) 
              
c Computing the rotation matrix

              beta = x(ilugan+1)
              gama = x(ilugan+2)
              teta = x(ilugan+3)
              
              call eulerrmat(beta,gama,teta,v1,v2,v3)  

              call comprindex(itype,imol,ifpen)

c Looping over the atoms of this molecule
  
              idatom = idfirst(itype) - 1
              do iatom = 1, natoms(itype) 
                 
                 call compicart(itype,imol,iatom,icart)
                 idatom = idatom + 1

c Computing the cartesian coordinates for this atom

                 call compcart(icart,xcart,
     +                    xbar,ybar,zbar,
     +                    coor(idatom,1),coor(idatom,2),coor(idatom,3),
     +                    v1,v2,v3)

c Adding to f the value relative to constraints for this atom

                 call comprest(xcart,restpars,
     +                    scale,scale2,
     +                    nratom,ityperest,iratom,
     +                    icart,fplus)

                 fpen(ifpen) = fpen(ifpen) + fplus
                 fresttmp = dmax1(fresttmp,fplus)
                 if(move) fatom(icart) = fatom(icart) + fplus

              end do 
           end do
!$OMP END DO


!$OMP ATOMIC
           frest = dmax1(frest,fresttmp)

        end if
      end do

!$OMP END PARALLEL
c      call system_clock(t2)
c      ftime = ftime + real(t2-t1)/rate
      
c Putting atoms in their boxes
      do itype = 1, ntype       
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

c Summing the penalities

      ftemp = 0.d0
      do itype = 1, ntotmol
               ftemp = ftemp + fpen(itype)
      end do
      
      f = f + ftemp


c ---

c Minimum distance function evaluation

      fdist = 1.d20
      tot_boxes = nboxes(1)*nboxes(2)*nboxes(3)

      if (tot_boxes .gt. nbp**3) write (*,*) "ERROR"

      if(.not. init1) then

c      call system_clock(t1,rate)

!$OMP PARALLEL PRIVATE(my_rank,num_proc,min_box,max_box,i,j,k,icart,
!$OMP+fdtmp,box_counter,ktemp,nl,itemp,ftemp,flast) DEFAULT(SHARED)

      ftemp = 0.d0
      fdtmp = fdist

!$OMP DO SCHEDULE(GUIDED)
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

c           Vector that keeps the value for this atom

               if(move) flast = ftemp


c           Interactions inside box

               ftemp = ftemp + fparc(icart,latomnext(icart),fdtmp)

c           Interactions of boxes that share faces

               ftemp = ftemp + fparc(icart,latomfirst(i+1,j,k),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i,j+1,k),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i,j,k+1),fdtmp)

c           Interactions of boxes that share axes

               ftemp = ftemp + fparc(icart,latomfirst(i+1,j+1,k),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i+1,j,k+1),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i+1,j-1,k),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i+1,j,k-1),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i,j+1,k+1),fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i,j+1,k-1),fdtmp)

c           Interactions of boxes that share vertices

               ftemp = ftemp + fparc(icart,latomfirst(i+1,j+1,k+1),
     +                               fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i+1,j+1,k-1),
     +                               fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i+1,j-1,k+1),
     +                               fdtmp)
               ftemp = ftemp + fparc(icart,latomfirst(i+1,j-1,k-1),
     +                               fdtmp)

c           If going to move bad molecules, update fatom

               if(move) fatom(icart) = fatom(icart) + f - flast

            end if

            icart = latomnext(icart)
         end do

         farray(box_counter) = ftemp
         ftemp = 0
      end do
!$OMP END DO

!$OMP ATOMIC
      fdist = dmin1(fdist,fdtmp)


!$OMP END PARALLEL

c      call system_clock(t2)
c      ftime = ftime + real(t2-t1)/rate

      ftemp = 0.d0
      do box_counter=1,tot_boxes
         ftemp = ftemp + farray(box_counter)
      end do

      f = f + ftemp

      end if

      return

      end

c
c Function that computes the main part of the objective function
c
      double precision function fparc(icart,firstjcart,fdtmp)

c Dimensions
      include 'sizes.i'

C     SCALAR ARGUMENTS
      integer icart,firstjcart

C     OUTPUT
      double precision fdtmp

C     LOCAL SCALARS
      integer jcart
      double precision a1,a2,a3,datom

C     COMMON BLOCK
      include 'molpa.i'

      fparc = 0.0d0
      jcart = firstjcart
      do while ( jcart .ne. 0 )
        if(comptype(ibtype(jcart))) then
          if(ibmol(icart).ne.ibmol(jcart).or.
     +       ibtype(icart).ne.ibtype(jcart)) then
            a1 = xcart(icart, 1)-xcart(jcart, 1) 
            a2 = xcart(icart, 2)-xcart(jcart, 2) 
            a3 = xcart(icart, 3)-xcart(jcart, 3) 
            datom = a1 * a1 + a2 * a2 + a3 * a3
            a1 = dmin1(datom - dism2, 0.d0)
            fparc = fparc + a1 * a1
            fdtmp = dmin1(datom,fdtmp)
          end if
        end if
        jcart = latomnext(jcart)
      end do

      end
