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

      double precision x(nn), f

      call feasyseq(x,f)

      return

      end

