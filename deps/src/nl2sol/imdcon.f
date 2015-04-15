      integer function imdcon(k)
c
      integer k
c
c  ***  return integer machine-dependent constants  ***
c
c     ***  k = 1 means return standard output unit number.   ***
c     ***  k = 2 means return alternate output unit number.  ***
c     ***  k = 3 means return  input unit number.            ***
c          (note -- k = 2, 3 are used only by test programs.)
c
      integer mdcon(3)
c      data mdcon(1)/6/, mdcon(2)/8/, mdcon(3)/5/
      data mdcon(1)/6/, mdcon(2)/6/, mdcon(3)/5/
c
      imdcon = mdcon(k)
      return
c  ***  last card of imdcon follows  ***
      end
