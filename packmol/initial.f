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
c Subroutine initial: Subroutine that reset parameters and
c                     builds the initial point
c

      subroutine initial(isem,randini,x,n,ntfix,fix,
     +                   moldy,chkgrad,nloop,discale,precision,sidemax,
     +                   movefrac,check)

      implicit none
      include 'sizes.i'
      include 'molpa.i'
      include 'gencan.i'

      integer isem, n, i, j, k, idatom, iatom, ilubar, icart, itype, 
     +        imol, ntry, ntfix, nb, iboxx, iboxy, iboxz, ifatom, 
     +        idfatom, iftype, jatom, nloop

      double precision x(nn), cmx, cmy, 
     +                 cmz, fx, xlength, dbox, rnd, discale, precision, 
     +                 movefrac, twopi, sidemax
      parameter(twopi = 2.d0 * 3.1415925655d0)
     
      logical fix, randini, moldy, chkgrad, hasfixed(nbp,nbp,nbp),
     +        overlap, movebadprint, hasbad, check

c We need to initialize the move logical variable

      move = .false.

c Default status of the function evaluation

      init1 = .false.

c Initialize the comptype logical array

      do i = 1, ntfix
        comptype(i) = .true.
      end do

c Penalty factors for the objective function relative to restrictions
c Default values: scale = 1.d2, scale2 = 1.d1

      scale = 1.d0
      scale2 = 1.d-2

c Move molecules to their center of mass (not for moldy)                                                                                   
      if(.not.moldy) call tobar(coor,ntype,natoms,idfirst)

c Compute maximum internal distance within each type of molecule

      do itype = 1, ntype
        dmax(itype) = 0.d0
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype) - 1
          do jatom = iatom + 1, natoms(itype)
            dmax(itype) = dmax1 ( dmax(itype),
     +             (coor(idatom+iatom,1)-coor(idatom+jatom,1))**2+
     +             (coor(idatom+iatom,2)-coor(idatom+jatom,2))**2+
     +             (coor(idatom+iatom,3)-coor(idatom+jatom,3))**2 )
          end do
        end do
        dmax(itype) = dsqrt(dmax(itype))
        write(*,*) ' Maximum internal distance of type ',itype,': ',
     +             dmax(itype)
        if(dmax(itype).eq.0.) dmax(itype) = 1.d0
      end do

c Maximum size of the system: if you system is very large (about
c 80 nm wide), increase the sidemax parameter.
c Otherwise, the packing can be slow and unsucesful

      cmxmin(1) = -sidemax
      cmymin(1) = -sidemax
      cmzmin(1) = -sidemax
      cmxmax(1) = sidemax
      cmymax(1) = sidemax
      cmzmax(1) = sidemax
      do i = 1, 3
        x(i) = 0.d0
        x(i+ntotmol*3) = 0.d0
      end do
      call restmol(1,0,n,x,fx,.true.,movefrac,precision,isem)
      sizemin(1) = x(1) - sidemax 
      sizemax(1) = x(1) + sidemax
      sizemin(2) = x(2) - sidemax
      sizemax(2) = x(2) + sidemax
      sizemin(3) = x(3) - sidemax
      sizemax(3) = x(3) + sidemax
      write(*,*) ' All atoms must be within these coordinates: '
      write(*,*) '  x: [ ', sizemin(1),', ', sizemax(1), ' ] '
      write(*,*) '  y: [ ', sizemin(2),', ', sizemax(2), ' ] '
      write(*,*) '  z: [ ', sizemin(3),', ', sizemax(3), ' ] '
      write(*,*) ' If the system is larger than this, increase the',
     +           ' sidemax parameter. '

c Create first aleatory guess

      do i = 1, n/2, 3
        x(i) = sizemin(1) + rnd(isem)*(sizemax(1)-sizemin(1))
        x(i+1) = sizemin(2) + rnd(isem)*(sizemax(2)-sizemin(2))
        x(i+2) = sizemin(3) + rnd(isem)*(sizemax(3)-sizemin(3))
      end do
      do i = n/2 + 1, n
        x(i) = twopi * rnd(isem)
      end do

c Compare analytical and finite-difference gradients

      if(chkgrad) then
        dbox = discale * dism + 0.01d0 * dism 
        do i = 1, 3
          xlength = sizemax(i) - sizemin(i)
          nb = int(xlength/dbox + 1.d0)  
          if(nb.gt.nbp) nb = nbp
          boxl(i) = dmax1(xlength/dfloat(nb),dbox)
          nboxes(i) = nb
        end do
        call compgrad(n,x)
        stop
      end if

c Performing some steps of optimization for the restrictions only

      write(*,*)
      write(*,905)
905   format( /, 13('-'),
     +       ' Building initial approximation ... ', 13('-'),/,/,
     +       '  Adjusting initial point to fit the constraints ')   
      init1 = .true.
      i = 0
      fx = 1.d0
700   format('  Packing:|0 ',tr39,'  10|')   
      hasbad = .true.
      do while( fx.gt.precision .and. i.le. (nloop/10-1) .and. hasbad)
        i = i + 1 
        write(*,700)
        call pgencan(n,x,fx)
        call feasy(x,fx)
        if(fx.gt.precision) then 
          write(*,701)'  Moving worst molecules ... ', i,' of ',nloop/10
701       format(a,i6,a,i6)
          movebadprint = .false.
          call movebad(n,x,fx,movefrac,precision,isem,
     +                 hasbad,movebadprint)
        end if
      end do
      write(*,*) 
      write(*,*) ' Restraint-only function value: ', fx
      init1 = .false.

      if(hasbad .and. fx.gt.precision) then
        write(*,*) ' ERROR: Packmol was unable to put the molecules'
        write(*,*) '        in the desired regions even without'
        write(*,*) '        considering distance tolerances. '
        write(*,*) '        Probably there is something wrong with'
        write(*,*) '        the constraints, since it seems that'
        write(*,*) '        the molecules cannot satisfy them at'
        write(*,*) '        at all. '
        write(*,*) '        Please check the spatial constraints and' 
        write(*,*) '        try again.'
        if ( i .ge. nloop/10-1 ) then
        end if
          write(*,*) ' >The maximum number of cycles (',nloop,
     +               ') was achieved.' 
          write(*,*) '  You may try increasing it with the',
     +               ' nloop keyword, as in: nloop 1000 '
        stop
      end if

c Rescaling sizemin and sizemax in order to build the patch of boxes

      write(*,*) ' Rescaling maximum and minimum coordinates... '
      do i = 1, 3
        sizemin(i) = 1.d20
        sizemax(i) = -1.d20
      end do                       
      
      icart = 0
      do itype = 1, ntfix
        do imol = 1, nmols(itype)
          do iatom = 1, natoms(itype) 
            icart = icart + 1
            sizemin(1) = dmin1(sizemin(1),xcart(icart,1))
            sizemin(2) = dmin1(sizemin(2),xcart(icart,2))
            sizemin(3) = dmin1(sizemin(3),xcart(icart,3))
            sizemax(1) = dmax1(sizemax(1),xcart(icart,1))
            sizemax(2) = dmax1(sizemax(2),xcart(icart,2))
            sizemax(3) = dmax1(sizemax(3),xcart(icart,3))
          end do 
        end do
      end do             

c Computing the size of the patches

      write(*,*) ' Computing size of patches... '
      dbox = discale * dism + 0.01d0 * dism 
      do i = 1, 3
        xlength = sizemax(i) - sizemin(i)
        nb = int(xlength/dbox + 1.d0)  
        if(nb.gt.nbp) nb = nbp
        boxl(i) = dmax1(xlength/dfloat(nb),dbox)
        nboxes(i) = nb
      end do

c Reseting latomfix array

      do i = 1, nbp
        do j = 1, nbp
          do k = 1, nbp
            latomfix(i,j,k) = 0
            hasfixed(i,j,k) = .false.
          end do
        end do
      end do   
 
c If there are fixed molecules, add them permanently to the latomfix array

      write(*,*) ' Add fixed molecules to permanent arrays... '
      if(fix) then
        icart = 0
        do iftype = ntype + 1, ntfix
          idfatom = idfirst(iftype) - 1
          do ifatom = 1, natoms(iftype)
            idfatom = idfatom + 1
            icart = icart + 1
            xcart(icart,1) = coor(idfatom,1)
            xcart(icart,2) = coor(idfatom,2)
            xcart(icart,3) = coor(idfatom,3)
            call setibox(xcart(icart,1),
     +                   xcart(icart,2),
     +                   xcart(icart,3),
     +                   sizemin,boxl,nboxes,iboxx,iboxy,iboxz)
            latomnext(icart) = latomfix(iboxx,iboxy,iboxz)
            latomfix(iboxx,iboxy,iboxz) = icart
            ibtype(icart) = iftype
            ibmol(icart) = 1
            hasfixed(iboxx,  iboxy,  iboxz  ) = .true.
            hasfixed(iboxx+1,iboxy,  iboxz  ) = .true.
            hasfixed(iboxx,  iboxy+1,iboxz  ) = .true.
            hasfixed(iboxx,  iboxy,  iboxz+1) = .true.
            hasfixed(iboxx+1,iboxy+1,iboxz  ) = .true.
            hasfixed(iboxx+1,iboxy,  iboxz+1) = .true.
            hasfixed(iboxx+1,iboxy-1,iboxz  ) = .true.
            hasfixed(iboxx+1,iboxy,  iboxz-1) = .true.
            hasfixed(iboxx,  iboxy+1,iboxz+1) = .true.
            hasfixed(iboxx,  iboxy+1,iboxz-1) = .true.
            hasfixed(iboxx+1,iboxy+1,iboxz+1) = .true.
            hasfixed(iboxx+1,iboxy+1,iboxz-1) = .true.
            hasfixed(iboxx+1,iboxy-1,iboxz+1) = .true.
            hasfixed(iboxx+1,iboxy-1,iboxz-1) = .true.
          end do
        end do
      end if

c Reseting mass centers to be within the regions

      write(*,*) ' Reseting center of mass... '
      do itype = 1, ntype
        cmxmin(itype) = 1.d20
        cmymin(itype) = 1.d20
        cmzmin(itype) = 1.d20
        cmxmax(itype) = -1.d20
        cmymax(itype) = -1.d20
        cmzmax(itype) = -1.d20
      end do

      icart = natfix
      do itype = 1, ntype
        do imol = 1, nmols(itype)
          cmx = 0.d0
          cmy = 0.d0
          cmz = 0.d0
          do iatom = 1, natoms(itype)
            icart = icart + 1
            cmx = cmx + xcart(icart,1)
            cmy = cmy + xcart(icart,2)
            cmz = cmz + xcart(icart,3)
          end do
          cmx = cmx / dfloat(natoms(itype))
          cmy = cmy / dfloat(natoms(itype))
          cmz = cmz / dfloat(natoms(itype))
          cmxmin(itype) = dmin1(cmxmin(itype),cmx)
          cmymin(itype) = dmin1(cmymin(itype),cmy)
          cmzmin(itype) = dmin1(cmzmin(itype),cmz)
          cmxmax(itype) = dmax1(cmxmax(itype),cmx)
          cmymax(itype) = dmax1(cmymax(itype),cmy)
          cmzmax(itype) = dmax1(cmzmax(itype),cmz)
        end do
      end do

c Building random initial point 

      write(*,*) ' Building random initial point ... '
      do i = n/2 + 1, n
        x(i) = twopi * rnd(isem)
      end do
      ilubar = 0
      do itype = 1, ntype
        do imol = 1, nmols(itype)
          fx = 1.d0
          ntry = 0
          overlap = .false.
          do while((overlap.or.fx.gt.precision).and.ntry.le.20) 
            ntry = ntry + 1
            x(ilubar+1) = cmxmin(itype) +
     +                    rnd(isem)*(cmxmax(itype)-cmxmin(itype))
            x(ilubar+2) = cmymin(itype) +
     +                    rnd(isem)*(cmymax(itype)-cmymin(itype))
            x(ilubar+3) = cmzmin(itype) +
     +                    rnd(isem)*(cmzmax(itype)-cmzmin(itype))
            if(fix) then
              call setibox(x(ilubar+1),
     +                     x(ilubar+2),
     +                     x(ilubar+3),
     +                     sizemin,boxl,nboxes,iboxx,iboxy,iboxz)
              if(hasfixed(iboxx,  iboxy,  iboxz  ).or.
     +           hasfixed(iboxx+1,iboxy,  iboxz  ).or.
     +           hasfixed(iboxx,  iboxy+1,iboxz  ).or.
     +           hasfixed(iboxx,  iboxy,  iboxz+1).or.
     +           hasfixed(iboxx+1,iboxy+1,iboxz  ).or.
     +           hasfixed(iboxx+1,iboxy,  iboxz+1).or.
     +           hasfixed(iboxx+1,iboxy-1,iboxz  ).or.
     +           hasfixed(iboxx+1,iboxy,  iboxz-1).or.
     +           hasfixed(iboxx,  iboxy+1,iboxz+1).or.
     +           hasfixed(iboxx,  iboxy+1,iboxz-1).or.
     +           hasfixed(iboxx+1,iboxy+1,iboxz+1).or.
     +           hasfixed(iboxx+1,iboxy+1,iboxz-1).or.
     +           hasfixed(iboxx+1,iboxy-1,iboxz+1).or.
     +           hasfixed(iboxx+1,iboxy-1,iboxz-1)) then
                overlap = .true.
              else
                overlap = .false.
              end if
            end if  
            if(.not.overlap) call restmol(itype,ilubar,n,x,fx,.false.,
     +                                    movefrac,precision,isem)
          end do
          ilubar = ilubar + 3
        end do
      end do

c Return with current random point (not default)

      if(randini) return
 
c Adjusting current point to fit the constraints

      init1 = .true.
      i = 0
      fx = 1.d0
      hasbad = .true.
      do while( fx.gt.precision .and. i.le. (nloop/10-1) .and. hasbad)
        i = i + 1 
        write(*,700)
        call pgencan(n,x,fx)
        call feasy(x,fx)
        if(fx.gt.precision) then
          write(*,701)'  Moving worst molecules ... ', i,' of ',nloop/10
          movebadprint = .false.
          call movebad(n,x,fx,movefrac,precision,isem,
     +                 hasbad,movebadprint)
        end if
      end do
      write(*,*) 
      write(*,*) ' Restraint-only function value: ', fx
      init1 = .false.

      write(*,907)
907   format(/,62('#'),/)

      return
      end

c
c Subroutine resetboxes
c

      subroutine resetboxes(nboxes,latomfirst,latomfix)
      
      implicit none
      include 'sizes.i'

      integer i,j,k
      integer nboxes(3)
      integer latomfirst(0:nbp+1,0:nbp+1,0:nbp+1),
     +        latomfix(nbp,nbp,nbp) 

c Reset boxes

      do i = 1, nboxes(1)
        do j = 1, nboxes(2)
          do k = 1, nboxes(3)
            latomfirst(i,j,k) = latomfix(i,j,k)
          end do
        end do
      end do

c Reset margins
      
      do j = 0, nboxes(2)+1
        do k = 0, nboxes(3)+1
          latomfirst(0,j,k) = 0
          latomfirst(nboxes(1)+1,j,k) = 0
        end do
      end do

      do i = 0, nboxes(1)+1
        do k = 0, nboxes(3)+1
          latomfirst(i,0,k) = 0
          latomfirst(i,nboxes(2)+1,k) = 0
        end do
      end do

      do i = 0, nboxes(1)+1
        do j = 0, nboxes(2)+1
          latomfirst(i,j,0) = 0
          latomfirst(i,j,nboxes(3)+1) = 0
        end do
      end do      

      return
      end


c 
c Random number generator
c 
      integer function mult( p, q) 
  
      Integer p, q, p0, p1, q0, q1 
  
      p1 = p/10000 
      p0 = mod(p,10000) 
      q1 = q/10000 
      q0 = mod(q,10000) 
      mult = mod( mod( p0*q1+p1*q0,10000)*10000+p0*q0,100000000) 
      return 
      end 
 
      double precision function rnd(sem) 
  
      integer sem, mult 
  
      sem = mod( mult( sem, 3141581) + 1, 100000000) 
      rnd = sem/100000000.0d0 
      return 
      end 

c
c subroutine tobar: moves molecules to their baricentres
c

      subroutine tobar(coor,ntype,natoms,idfirst)
      
      implicit none
      include 'sizes.i'
      integer ntype, natoms(maxtype), idfirst(maxtype), idatom,
     +        itype, iatom
      double precision xcm, ycm, zcm, coor(maxatom,3)

      do itype = 1, ntype
        idatom = idfirst(itype) - 1
        xcm = 0.d0
        ycm = 0.d0
        zcm = 0.d0
        do iatom = 1, natoms(itype)
          idatom = idatom + 1
          xcm = xcm + coor(idatom,1)
          ycm = ycm + coor(idatom,2)
          zcm = zcm + coor(idatom,3)
        end do
        xcm = xcm / natoms(itype)
        ycm = ycm / natoms(itype)
        zcm = zcm / natoms(itype)
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype)
          idatom = idatom + 1
          coor(idatom,1) = coor(idatom,1) - xcm
          coor(idatom,2) = coor(idatom,2) - ycm
          coor(idatom,3) = coor(idatom,3) - zcm
        end do
      end do

      return                                                 
      end  

c
c Subroutine setibox: set box index for given coordinates
c 

      subroutine setibox(x,y,z,sizemin,boxl,nboxes,iboxx,iboxy,iboxz)

      implicit none
      double precision x, y, z, sizemin(3), boxl(3), xtemp, ytemp, ztemp
      integer nboxes(3), iboxx, iboxy, iboxz

      xtemp = x - sizemin(1) 
      ytemp = y - sizemin(2)
      ztemp = z - sizemin(3)
      iboxx = int(xtemp/boxl(1)) + 1
      iboxy = int(ytemp/boxl(2)) + 1
      iboxz = int(ztemp/boxl(3)) + 1
      if(xtemp.le.0) iboxx = 1
      if(ytemp.le.0) iboxy = 1
      if(ztemp.le.0) iboxz = 1 
      if(iboxx.gt.nboxes(1)) iboxx = nboxes(1)
      if(iboxy.gt.nboxes(2)) iboxy = nboxes(2)
      if(iboxz.gt.nboxes(3)) iboxz = nboxes(3)

      return
      end

c
c subroutine restmol: either compute the restraint function
c                     value for a single molecule or solve
c                     the problem of puting this molecule
c                     in the restraint region
c  

      subroutine restmol(itype,ilubar,n,x,fx,solve,movefrac,
     +                   precision,isem)

      implicit none
      include 'sizes.i'
      include 'molpa.i'
      include 'gencan.i'

      integer n, nsafe, ntotsafe, itype, i, ilubar, nmoltype, isem,
     +        ip1, ip2
      double precision xmol(nn), x(nn), fx, movefrac, precision
      logical solve, compsafe(maxtype), initsafe
c      logical movebadprint, hasbad 

c Saving global problem variables

      nsafe = n
      ntotsafe = ntotmol
      nmoltype = nmols(itype)
      do i = 1, ntype
        compsafe(i) = comptype(i)
      end do
      initsafe = init1

c Preparing system to solve for this molecule

      n = 6
      ntotmol = 1      
      nmols(itype) = 1
      xmol(1) = x(ilubar+1)
      xmol(2) = x(ilubar+2)
      xmol(3) = x(ilubar+3)
      xmol(4) = x(ilubar+ntotsafe*3+1)
      xmol(5) = x(ilubar+ntotsafe*3+2)
      xmol(6) = x(ilubar+ntotsafe*3+3)
      do i = 1, ntype
        if(i.eq.itype) then
          comptype(i) = .true.
        else
          comptype(i) = .false.
        end if
      end do
      init1 = .true.
      
c If not going to solve the problem, compute energy and return

      if(.not.solve) then
        call feasy(xmol,fx)

c Otherwise, put this molecule in its constraints

      else
        ip1 = iprint1
        ip2 = iprint2
        iprint1 = 0
        iprint2 = 0
        call pgencan(n,xmol,fx)
c        movebadprint = .false.
c        call movebad(n,xmol,fx,movefrac,precision,isem,
c     +                hasbad,movebadprint)
        iprint1 = ip1
        iprint2 = ip2
      end if       

c Restoring original problem data

      ntotmol = ntotsafe
      n = nsafe
      nmols(itype) = nmoltype
      x(ilubar+1) = xmol(1) 
      x(ilubar+2) = xmol(2) 
      x(ilubar+3) = xmol(3) 
      x(ilubar+ntotmol*3+1) = xmol(4) 
      x(ilubar+ntotmol*3+2) = xmol(5) 
      x(ilubar+ntotmol*3+3) = xmol(6) 
      do i = 1, ntype
        comptype(i) = compsafe(i)
      end do
      init1 = initsafe

      return
      end



