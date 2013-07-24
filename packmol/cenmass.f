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
c Subroutine cenmass
c
c            Computes the center of mass of free molecules and
c            for fixed molecules, if required. 
c

      subroutine cenmass(coor,amass,
     +                   ntype,nlines,
     +                   idfirst,natoms,
     +                   keyword,linestrut)

      implicit none
      include 'sizes.i'

      character*200 keyword(maxlines,maxkeywords)

      double precision cm(maxtype,3), totm(maxtype)
      double precision amass(maxatom)
      double precision coor(maxatom,3)

      integer linestrut(maxtype,2)  
      integer k, iline, nlines
      integer itype, ntype, iatom, idatom
      integer idfirst(maxtype), natoms(maxtype)

      logical domass(maxtype)

c Setting the molecules for which the center of mass is computed

      do itype = 1, ntype
        domass(itype) = .true.
      end do

      do iline = 1, nlines
        if(keyword(iline,1).eq.'fixed') then
          do itype = 1, ntype
            if(iline.gt.linestrut(itype,1).and.
     +         iline.lt.linestrut(itype,2)) then
              domass(itype) = .false.
            end if
          end do
        end if
      end do
            
      do iline = 1, nlines
        if(keyword(iline,1).eq.'centerofmass'.or.
     +     keyword(iline,1).eq.'center') then
          do itype = 1, ntype
            if(iline.gt.linestrut(itype,1).and.
     +         iline.lt.linestrut(itype,2)) then
              domass(itype) = .true.
            end if
          end do
        end if
      end do

c Computing the center of mass

      do itype = 1, ntype 
        do k = 1, 3 
          cm(itype, k) = 0.d0 
        end do 
      end do 
 
      do itype = 1, ntype 
        totm(itype) = 0.d0 
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype) 
          idatom = idatom + 1
          totm(itype) = totm(itype) + amass(idatom) 
        end do 
      end do 
 
      do itype = 1, ntype 
        idatom = idfirst(itype) - 1
        do iatom = 1, natoms(itype)
          idatom = idatom + 1 
          do k = 1, 3 
            cm(itype, k) = cm(itype, k)  
     +                     + coor(idatom, k)*amass(idatom) 
          end do 
        end do 
        do k = 1, 3 
          cm(itype, k) = cm(itype, k) / totm(itype) 
        end do 
      end do  


c Putting molecules in their center of mass

      do itype = 1, ntype
        if(domass(itype)) then
          idatom = idfirst(itype) - 1
          do iatom = 1, natoms(itype)
            idatom = idatom + 1
            do k = 1, 3
              coor(idatom, k) = coor(idatom, k) - cm(itype, k)
            end do
          end do
        end if
      end do

      return
      end
