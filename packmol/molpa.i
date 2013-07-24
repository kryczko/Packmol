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
      double precision xcart(maxatom,3) 
      double precision coor(maxatom,3) 
      double precision restpars(maxrest,9) 
      double precision sizemin(3),sizemax(3)
      double precision boxl(3)
      double precision fdist, frest 
      double precision rot_bound(maxtype,3,2)
      
      double precision dism, dism2
      double precision scale, scale2
      double precision fatom(maxatom)
      double precision dmax(maxtype)
      double precision cmxmin(maxtype),cmymin(maxtype),cmzmin(maxtype)
      double precision cmxmax(maxtype),cmymax(maxtype),cmzmax(maxtype)

      integer nmols(maxtype)    
      integer natoms(maxtype) 
      integer idfirst(maxtype)
      integer nratom(maxatom)   
      integer iratom(maxatom,mrperatom) 
      integer ityperest(maxrest)  
      integer nboxes(3)  
      integer ibmol(maxatom)  
      integer ibtype(maxatom)  

      integer ntotmol, ntype, natfix, ntotat

      logical constrain_rot(maxtype,3)
      logical comptype(maxtype)
      logical init1, move

      integer latomnext(maxatom),latomfirst(0:nbp+1,0:nbp+1,0:nbp+1),
     +        latomfix(nbp,nbp,nbp)

      common/evalfg/xcart,fatom,coor,restpars,sizemin,sizemax,boxl,
     +              cmxmin, cmymin, cmzmin, cmxmax, cmymax, cmzmax,
     +              dmax,fdist,frest,rot_bound,
     +              dism,dism2,scale,scale2,
     +              nmols,natoms,idfirst,nratom,iratom,ityperest,
     +              latomnext,latomfirst,latomfix,
     +              nboxes,ibmol,ibtype,
     +              ntotat,ntotmol,ntype,natfix,
     +              constrain_rot,comptype,init1,move

