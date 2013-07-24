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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
c
c Packmol: A package for building initial configurations for
c molecular dynamics simulations, to be published, 2008.
c
c http://www.ime.unicamp.br/~martinez/packmol
c
c Usage (see the page above for further information):
c
c ./packmol < inputfile.inp
c
c References:
c
c L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez,
c PACKMOL: A package for building initial configurations for
c molecular dynamics simulations, J. Comp. Chem. 30:2157-2164, 2009.
c
c J. M. Martinez and L. Martinez, 
c Packing optimization for the automated generation of complex
c system's initial configurations for molcular dynamics and
c docking. J. Comp. Chem. 24:819-825, 2003.
c
c This version of Packmol uses the optimization method GENCAN which
c is a part of the TANGO (Trustable Algorithms for Nonlinear General
c Optimization) project.
c Reference:
c E. G. Birgin, J. M. Martinez, Comp. Opt. Appl. 23:101-125, 2002.
c http://www.ime.usp.br/~egbirgin/tango
c
c

      program packmol

      implicit none
      include 'sizes.i'
      include 'molpa.i'
      include 'gencan.i'

      integer irestline(maxrest)
      integer linestrut(maxtype,2)
      integer itype, nrest, irest, idatom, iatom
      integer ntemp, ntfix, idtemp, nmtemp, natemp, nlines
      integer linesttmp1, linesttmp2, jtype
      integer ntmol, n, iftype, icart, imol 
      integer i, iline, iiatom, iat, iirest, iratcount
      integer isem
      integer nloop, loop
      integer maxcon(maxatom)
      integer ntcon(9), nconnect(maxatom,8) 
      integer ntottemp, ilubar, ilugan  
      integer resnumbers(maxtype), resntemp
      integer charl, writeout
      integer input_itype(maxtype)
c      integer iargc, narg, charl, writeout
      
      double precision v1(3),v2(3),v3(3)
      double precision x(nn), xfull(nn), xbest(nn)
      double precision disini
      double precision cmx, cmy, cmz, beta, gama, teta
      double precision xtemp, ytemp, ztemp
      double precision fx, bestf, flast, fout
      double precision fimp, fimprov, precision
      double precision amass(maxatom), charge(maxatom)
      double precision discale, movefrac
      double precision add_sides_fix
      double precision sidemax
      double precision pi
      parameter(pi=3.141592653589793d0)

      real etime, tarray(2), time0
      
      character*200 keyword(maxlines,maxkeywords)
      character*200 record
      character*200 name(maxtype)
      character*80 pdbfile(maxtype), xyzfile
      character*3 ele(maxatom)
      character*200 xyzout        

      logical fix,fixed(maxtype),fixtmp,randini,check,chkgrad
      logical pdb,tinker,xyz,moldy,rests,writebad
      logical add_amber_ter, add_box_sides
      logical movebadprint, hasbad
      logical thisisfixed(maxtype)

c Start time computation

      time0 = etime(tarray)

c If the user tried to run without redirection, output error
c
c      narg = iargc()
c      call getarg(1,record)
c      if(record.gt." ") then
c        if(narg.gt.0) then
c          if(record.le."0") then
c            record = "inputfile.inp"
c          end if
c          write(*,*)
c          write(*,*) ' Packmol must be run with input file redirection:'
c          write(*,*) ' packmol < ',record(1:charl(record))
c          write(*,*)
c          stop
c        end if
c      end if
c
c Printing title

1     format( 62('#'), /,/
     +        ' PACKMOL - Packing optimization for the automated', /
     +        ' generation of starting configurations for',        /
     +        ' molecular dynamics. ',/
     +        ' ',/
     +        t42,' Version 13.112 ',/
     +        ,/,62('#'),               /,/)
      write(*,1)

c Reading input file

      call getinp(dism,precision,sidemax,
     +            ntype,nlines,nrest,
     +            natoms,idfirst,nconnect,maxcon,nmols,
     +            isem,
     +            discale,nloop,
     +            irestline,ityperest,linestrut,
     +            coor,amass,charge,restpars,
     +            pdbfile,name,ele,keyword,
     +            xyzout,writeout,writebad,
     +            tinker,pdb,xyz,moldy,check,chkgrad,
     +            randini,resnumbers,movefrac,
     +            add_amber_ter,add_box_sides,add_sides_fix)

c Put molecules in their center of mass

      call cenmass(coor,amass,
     +             ntype,nlines,
     +             idfirst,natoms,
     +             keyword,linestrut)
 
c Computing the total number of atoms
     
      ntotat = 0
      do itype = 1, ntype
        ntotat = ntotat + natoms(itype) * nmols(itype)
      end do              
      write(*,*) ' Total number of atoms: ', ntotat

      if(ntotat.gt.maxatom) then
        write(*,*)' ERROR: Total number of atoms greater than maxatom.'
        write(*,*)'        Change the maxatom (sizes.i file) '
        stop
      end if

c Put fixed molecules in the specified position

      do itype = 1, ntype
        fixed(itype) = .false.
      end do

      do irest = 1, nrest
        if(ityperest(irest).eq.1) then
          do itype = 1, ntype
            if(irestline(irest).gt.linestrut(itype,1).and.
     +         irestline(irest).lt.linestrut(itype,2)) then
              cmx = restpars(irest,1) 
              cmy = restpars(irest,2)
              cmz = restpars(irest,3)    
              beta = restpars(irest,4) 
              gama = restpars(irest,5) 
              teta = restpars(irest,6) 

c Compute rotation matrix from euler angles

              call eulerfixed(beta,gama,teta,v1,v2,v3)                 

              idatom = idfirst(itype) - 1
              do iatom = 1, natoms(itype)
                idatom = idatom + 1
                xtemp =   coor(idatom,1)*v1(1) 
     +                  + coor(idatom,2)*v2(1) 
     +                  + coor(idatom,3)*v3(1) 
                ytemp =   coor(idatom,1)*v1(2) 
     +                  + coor(idatom,2)*v2(2) 
     +                  + coor(idatom,3)*v3(2) 
                ztemp =   coor(idatom,1)*v1(3) 
     +                  + coor(idatom,2)*v2(3) 
     +                  + coor(idatom,3)*v3(3) 
                coor(idatom, 1) = xtemp + cmx
                coor(idatom, 2) = ytemp + cmy
                coor(idatom, 3) = ztemp + cmz 
              end do
              record = name(itype)
              write(*,*) ' Molecule ',record(1:charl(record)),
     +                   '(',itype,') will be fixed.' 
              fixed(itype) = .true.
              if(nmols(itype).gt.1) then
                write(*,*)' ERROR: You cannot set number > 1',
     +                    ' for fixed molecules. '
                stop
              end if
            end if
          end do
        end if
      end do 

c Reseting parameters for removing the fixed molecules

      fix = .false.
      ntemp = 0
      do itype = 1, ntype

c input_itype and thisisfixed vectors are used only to preserve the
c order of input in the output files

        input_itype(itype) = itype
        if(fixed(itype)) then
          fix = .true.
          thisisfixed(itype) = .true.
        else
          ntemp = ntemp + 1
          thisisfixed(itype) = .false.
        end if
      end do
      ntfix = ntype
      ntype = ntemp     

      do i = 1, ntfix - ntype 
        do itype = 1, ntfix - 1
          if(fixed(itype)) then
            record = name(itype)
            fixtmp = fixed(itype)
            idtemp = idfirst(itype)
            nmtemp = nmols(itype)
            natemp = natoms(itype)
            resntemp = resnumbers(itype)
            if(pdb) xyzfile = pdbfile(itype)
            linesttmp1 = linestrut(itype,1)
            linesttmp2 = linestrut(itype,2)
            jtype = itype + 1
            if(.not.fixed(jtype)) then
              name(itype) = name(jtype)
              name(jtype) = record(1:10)
              idfirst(itype) = idfirst(jtype)
              idfirst(jtype) = idtemp
              fixed(itype) = fixed(jtype)
              fixed(jtype) = fixtmp
              nmols(itype) = nmols(jtype)
              nmols(jtype) = nmtemp
              natoms(itype) = natoms(jtype)
              natoms(jtype) = natemp
              resnumbers(itype) = resnumbers(jtype)
              resnumbers(jtype) = resntemp
              if(pdb) then
                pdbfile(itype) = pdbfile(jtype) 
                pdbfile(jtype) = xyzfile
              end if
              linestrut(itype,1) = linestrut(jtype,1)
              linestrut(itype,2) = linestrut(jtype,2)
              linestrut(jtype,1) = linesttmp1
              linestrut(jtype,2) = linesttmp2
            end if
          end if
        end do
      end do
 
c Computing the number of variables
c
c ntype: 1...ntype (counter for the number of free structures)
c
c ntfix: 1...ntype...ntfix (counter for the total number of structures)
c
      ntmol = 0
      do itype = 1, ntfix
        ntmol = ntmol + nmols(itype)
      end do
      ntotmol = 0 
      do itype = 1, ntype 
        ntotmol = ntotmol + nmols(itype)       
      end do     
      n = ntotmol * 6
      write(*,*) ' Total number of molecules: ', ntmol
      write(*,*) ' Number of fixed molecules: ', ntmol - ntotmol
      write(*,*) ' Number of free molecules: ', ntotmol
      write(*,*) ' Number of variables: ', n 

c Computing the total number of fixed atoms

      natfix = 0
      if(fix) then
        do iftype = ntype + 1, ntfix
          natfix = natfix + natoms(iftype)
        end do
      end if       
      write(*,*) ' Total number of fixed atoms: ', natfix

c Setting the array that contains the restrictions per atom

      icart = natfix
      do itype = 1, ntype
        rests = .false.
        do imol = 1, nmols(itype)
          idatom = idfirst(itype) - 1      
          do iatom = 1, natoms(itype) 
            icart = icart + 1
            idatom = idatom + 1
            nratom(icart) = 0
            iratcount = 0
            do i = 1, mrperatom
              iratom(icart,i) = 0
            end do
            iline = linestrut(itype,1)
            do while(iline.lt.linestrut(itype,2))
              iline = iline + 1
              if(keyword(iline,1).eq.'atoms') then
                iiatom = -1
                do iat = 2, maxkeywords
                  read(keyword(iline,iat),*,err=130,end=130) iiatom
                  if(iatom.eq.iiatom) goto 130
                end do
130             continue
                do while(keyword(iline,1).ne.'end'.and.
     +                   keyword(iline,2).ne.'atoms')
                  iline = iline + 1
                  if(iatom.eq.iiatom) then
                    if(keyword(iline,1).eq.'inside'.or.
     +                 keyword(iline,1).eq.'outside'.or.
     +                 keyword(iline,1).eq.'over'.or.
     +                 keyword(iline,1).eq.'below') then
                      nratom(icart) = nratom(icart) + 1
                      iratcount = iratcount + 1
                      do irest = 1, nrest
                        if(irestline(irest).eq.iline) iirest = irest
                      end do
                      iratom(icart,iratcount) = iirest
                    end if
                  end if
                end do
                iline = iline - 1
              else if(keyword(iline,1).eq.'inside'.or.
     +                keyword(iline,1).eq.'outside'.or.
     +                keyword(iline,1).eq.'over'.or.
     +                keyword(iline,1).eq.'below') then
                nratom(icart) = nratom(icart) + 1    
                iratcount = iratcount + 1
                do irest = 1, nrest
                  if(irestline(irest).eq.iline) iirest = irest
                end do
                iratom(icart,iratcount) = iirest
              end if
            end do
            if(nratom(icart).gt.0) rests = .true.
          end do 
          if(.not.rests) then
            write(*,*) ' ERROR: Some molecule has no geometrical',
     +                 ' restriction defined: nothing to do.'
            stop
          end if
        end do
      end do

c Read the constraints to rotations about axis, if set

      do itype = 1, ntype
        constrain_rot(itype,1) = .false.
        constrain_rot(itype,2) = .false.
        constrain_rot(itype,3) = .false.
        iline = linestrut(itype,1)
        do while(iline.lt.linestrut(itype,2))
          iline = iline + 1
          if(keyword(iline,1).eq.'constrain_rotation') then
            if(iline.gt.linestrut(itype,1).and.
     +         iline.lt.linestrut(itype,2)) then

c Note that for movable molecules, teta is a rotation on the x-axis,
c                                  gama is a rotation on the z-axis,
c                                  beta is a rotation on the y-axis
c                                  (see eulerrmat routine)

              if(keyword(iline,2).eq.'x') then
                constrain_rot(itype,3) = .true.
                read(keyword(iline,3),*) rot_bound(itype,3,1)
                read(keyword(iline,4),*) rot_bound(itype,3,2)
                rot_bound(itype,3,1) = rot_bound(itype,3,1)*pi/180.d0
                rot_bound(itype,3,2) = rot_bound(itype,3,2)*pi/180.d0
      
                write(*,*) ' Rotations about x axis of molecules of ',
     +          ' type ', itype, ' will be constrained. '
              end if
              if(keyword(iline,2).eq.'y') then
                constrain_rot(itype,1) = .true.
                read(keyword(iline,3),*) rot_bound(itype,1,1)
                read(keyword(iline,4),*) rot_bound(itype,1,2)
                rot_bound(itype,1,1) = rot_bound(itype,1,1)*pi/180.d0
                rot_bound(itype,1,2) = rot_bound(itype,1,2)*pi/180.d0

                write(*,*) ' Rotations about y axis of molecules of ',
     +          ' type ', itype, ' will be constrained. '
              end if
              if(keyword(iline,2).eq.'z') then
                constrain_rot(itype,2) = .true.
                read(keyword(iline,3),*) rot_bound(itype,2,1)
                read(keyword(iline,4),*) rot_bound(itype,2,2)
                rot_bound(itype,2,1) = rot_bound(itype,2,1)*pi/180.d0
                rot_bound(itype,2,2) = rot_bound(itype,2,2)*pi/180.d0

                write(*,*) ' Rotations about z axis of molecules of ',
     +          ' type ', itype, ' will be constrained. '
              end if
            end if
          end if
        end do
      end do
 
c If there are no variables (only fixed molecules, stop)

      if(n.eq.0) then
        call output(x,amass,
     +              irestline,linestrut,maxcon,ntcon,nconnect,
     +              nrest,ntfix,resnumbers,
     +              ele,pdbfile,xyzout,name,
     +              pdb,tinker,xyz,moldy,fix,
     +              add_amber_ter,add_box_sides,add_sides_fix,
     +              input_itype,thisisfixed)
        write(*,908)
        write(*,*) ' There are only fixed molecules, therefore '
        write(*,*) ' there is nothing to do. '
        write(*,*) ' The output file contains the fixed molecule '
        write(*,*) ' in the desired position. '
        write(*,908)
        stop
      end if
  
c
c (Re)setting parameters and building initial point
c

      call initial(isem,randini,x,n,ntfix,fix,moldy,
     +             chkgrad,nloop,discale,precision,sidemax,
     +             movefrac,check)

c Computing the energy at the initial point

      disini = dism
      dism2 = dism * dism
      call feasy(x,fx)
      write(*,*) ' Objective function at initial point: ', fx
      bestf = fx
      flast = fx
      fout = fx
      do i = 1, n
        xbest(i) = x(i)
      end do

      if(check) then
        call output(x,amass,
     +              irestline,linestrut,maxcon,ntcon,nconnect,
     +              nrest,ntfix,resnumbers,
     +              ele,pdbfile,xyzout,name,
     +              pdb,tinker,xyz,moldy,fix,
     +              add_amber_ter,add_box_sides,add_sides_fix,
     +              input_itype,thisisfixed)
        write(*,*) ' Wrote initial point to output file: ',
     +             xyzout(1:charl(xyzout)) 
        stop
      end if

c
c Main loop: first pack types of molecules separately, then
c pack all molecules together
c

      do i = 1, nn
        xfull(i) = x(i)
      end do
      ntemp = n
      ntottemp = ntotmol
      if(ntype.eq.1) then
        itype = 1
      else 
        itype = 0
      end if

      do while(itype.le.ntype+1)
        itype = itype + 1
 
c Use slightly larger tolerance than required to improove convergence

        dism = discale * disini
        dism2 = dism * dism
       
        if(itype.le.ntype) then
          if(nmols(itype).eq.1) itype = itype + 1
        end if

c Adjusting parameters for packing only this type

        if(itype.le.ntype) then
          write(*,*)
          write(*,908)
          write(*,*)
          write(*,*) ' Packing molecules of type ', itype
          write(*,*)
          write(*,908)
908       format(62('#'))
          do i = 1, ntype
            if(i.eq.itype) then
              comptype(i) = .true.
            else
              comptype(i) = .false.
            end if
          end do
          n = nmols(itype) * 6
          ntotmol = nmols(itype)
          ilubar = 0
          do i = 1, itype - 1
            ilubar = ilubar + nmols(i) * 3
          end do
          ilubar = ilubar + 1
          ilugan = ntemp/2 + ilubar 
          do i = 1, n / 2
            x(i) = xfull(ilubar)
            x(i+n/2) = xfull(ilugan)
            ilubar = ilubar + 1
            ilugan = ilugan + 1
          end do
        end if

c If itype=ntype+1 restore original vectors and pack all molecules

        if(itype.eq.ntype+1) then
          n = ntemp 
          ntotmol = ntottemp
          do i = 1, n
            x(i) = xfull(i)
          end do
          do itype = 1, ntype
            comptype(itype) = .true.
          end do
          if(ntype.gt.1) then
            write(*,*)
            write(*,908)
            write(*,*)
            write(*,*)' Solving the problem for all molecules together.'
            write(*,*)
            write(*,908)
            write(*,*)
          end if
        end if
 
        loop = -1
        do while(loop.le.nloop)
        loop = loop + 1

c Reseting the parameters relative to the improvement of the function
           
        if(loop.eq.0) then
          fimp = 1.d99
          fimprov = fimp
          call feasy(x,fx)
          bestf = fx
          flast = fx
        end if

c Moving bad molecules

        if(fimp.le.10.d0) then
          movebadprint = .true.
          call movebad(n,x,fx,movefrac,precision,isem,
     +                 hasbad,movebadprint)
          flast = fx
        end if

        if(loop.eq.nloop.and.itype.eq.ntype+1) then
          write(*,*)' STOP: Maximum number of GENCAN loops achieved.'
          call checkpoint(n,xbest,amass,
     +                    nrest,ntfix,nloop,
     +                    irestline,linestrut,maxcon,ntcon,nconnect,
     +                    ele,pdbfile,xyzout,name,
     +                    pdb,tinker,xyz,moldy,fix,
     +                    movefrac,precision,isem,resnumbers,
     +                    add_amber_ter,add_box_sides,add_sides_fix,
     +                    input_itype,thisisfixed)
          stop
        end if

        write(*,905) loop, dism   
905     format( /, 17('-'),' Starting GENCAN loop(',i4,') ',17('-'),/
     +  '                                         Tolerance:',f10.2) 

C CALL GENCAN

        write(*,700)
700     format('  Packing:|0 ',tr39,'  10|')
        call pgencan(n,x,fx)

c Compute the statistics of the last optimization loop

        call feasy(x,fx)
        if(bestf.gt.0.d0) fimprov = -100.d0 * (fx - bestf) / bestf
        if(bestf.eq.0.d0) fimprov = 100.d0
        if(flast.gt.0.d0) fimp = -100.d0 * (fx - flast) / flast
        if(flast.eq.0.d0) fimp = 100.d0
        fimp = dmin1(99.99d0,dmax1(-99.99d0,fimp))
        fimprov = dmin1(99.99d0,dmax1(-99.99d0,fimprov))

        write(*,907) fx, bestf, fimprov, fimp, 
     +               dsqrt(fdist), frest
907     format(                                                       /
     *         '  Function value from last GENCAN loop: f = ', e10.5, /
     *         '  Best function value before: f = ', e10.5,           /
     *         '  Improvement from best function value: ', f8.2, ' %',/
     *         '  Improvement from last loop: ', f8.2, ' %',          /
     *         '  Minimum distance between atoms: ', f12.6,           /
     *         '  Maximum violation of the constraints: ', e10.5,     /
     *         62('-'),/)

        flast = fx

c If the distance between molecules is satisfactory, restore disini

        if(dsqrt(fdist).ge.disini) then
          dism = disini
          dism2 = dism * dism
        end if

c Updating best point

        if(fx.le.bestf) then
          bestf = fx 
          if(itype.eq.ntype+1) then
            do i = 1, n
              xbest(i) = x(i)
            end do
          end if
        end if

c Writing output file

        if(itype.eq.ntype+1) then  
          if( ( ( disini-dsqrt(fdist).lt.precision.and.
     +          frest.lt.precision ) .or. bestf.lt.precision )
     +      .or. ( mod(loop+1,writeout).eq.0 .and.
     +           ( writebad .or. bestf.lt.fout ) ) ) then
            fout = bestf
            if(writebad) then
              call output(x,amass,
     +                    irestline,linestrut,maxcon,ntcon,nconnect,
     +                    nrest,ntfix,resnumbers,
     +                    ele,pdbfile,xyzout,name,
     +                    pdb,tinker,xyz,moldy,fix,
     +                    add_amber_ter,add_box_sides,add_sides_fix,
     +                    input_itype,thisisfixed)
              write(*,*) ' Current point written to file: ', 
     +                   xyzout(1:charl(xyzout)) 
            else
              call output(xbest,amass,
     +                    irestline,linestrut,maxcon,ntcon,nconnect,
     +                    nrest,ntfix,resnumbers,
     +                    ele,pdbfile,xyzout,name,
     +                    pdb,tinker,xyz,moldy,fix,
     +                    add_amber_ter,add_box_sides,add_sides_fix,
     +                    input_itype,thisisfixed)
              write(*,*) ' Best solution written to file: ', 
     +                   xyzout(1:charl(xyzout))
            end if
          end if
        end if

c When the solution is found, print success information and stop

        if((disini-dsqrt(fdist).lt.precision.and.
     +      frest.lt.precision).or.
     +      bestf.lt.precision) then

          if(itype.le.ntype) then
            write(*,908)
            write(*,*)' Packing solved for molecules of type', itype
            write(*,*)' Objective function value: ', bestf
            write(*,*)' Minimum distance between atoms: ',dsqrt(fdist)
            write(*,*)' Max. constraint violation: ', frest
            write(*,908)
            loop = nloop + 1      
          else
909         format(/, 62('#'),/,                               /,
     +        t27, ' Success! ',                               /,
     +        t10, ' Final objective function value: ', e10.5, /,
     +        t10, ' Minimum distance between atoms: ', f10.6, /,
     +        t10, ' Maximum violation of the constraints: ', e10.5,/,
     +        62('-'), /,
     +        ' Please cite this work if Packmol was useful: ',/,
     +   ' L. Martinez, R. Andrade, E. G. Birgin, J. M. Martinez, ',/,
     +   ' PACKMOL: A package for building initial configurations ',/,
     +   ' for molecular dynamics simulations. ',/,
     +   ' Journal of Computational Chemistry, 30:2157-2164,2009.',
     +        /,/,62('#'),/)
            write(*,909) bestf, dsqrt(fdist), frest
            write(*,*) '  Running time: ', 
     +                    etime(tarray) - time0,' seconds. ' 
            stop 
          end if
        end if

c End do of loop do:
        end do

        if(itype.le.ntype) then
          ilubar = 0
          do i = 1, itype - 1
            ilubar = ilubar + nmols(i)*3
          end do
          ilubar = ilubar + 1
          ilugan = ntemp/2 + ilubar
          do i = 1, n/2
            xfull(ilubar) = x(i)
            xfull(ilugan) = x(i+n/2)
            ilubar = ilubar + 1
            ilugan = ilugan + 1
          end do
        end if

c End of itype do:
      end do

      end
