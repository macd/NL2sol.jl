      subroutine xinit(p, x, nex)
c
c     *****parameters...
c
      integer nex, p
      double precision x(p)
c
c     ..................................................................
c
c     *****purpose...
c     this routine initializes the solution vector x according to
c     the initial values for the various test functions given in
c     references (1), (2), and (3).
c     subroutines testr and testj.  (see testr for references.)
c
c     *****parameter description...
c     on input...
c
c        nex is the test problem number.
c
c        p is the number of parameters.
c
c     on output...
c
c        x is the initial guess to the solution.
c
c     *****application and usage restrictions...
c     this routine is called by nltest.
c
c     ..................................................................
c
c     *****local variables...
      integer i
      double precision pp1inv
c     *****intrinsic functions...
c/+
      real float
      double precision dble
c/
      dfloat(ii) = dble(float(ii))
c
      go to (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
     1   1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
     2   2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 1900, 2100,
     3   2500, 1300, 1400, 1500, 1600),nex
c
c  ***  rosenbrock  ***
 100  x(1) = -1.2d0
      x(2) = 1.0d0
      go to 9999
c  ***  helix  ***
 200  x(1) = -1.0d0
      x(2) = 0.0d0
      x(3) = 0.0d0
      go to 9999
c  *** singular  ***
 300  x(1) = 3.0d0
      x(2) = -1.0d0
      x(3) = 0.0d0
      x(4) = 1.0d0
      go to 9999
c  ***  woods  ***
 400  x(1) = -3.0d0
      x(2) = -1.0d0
      x(3) = -3.0d0
      x(4) = -1.0d0
      go to 9999
c  ***  zangwill  ***
 500  x(1) = 1.0d2
      x(2) = -1.0d0
      x(3) = 2.5d0
      go to 9999
c  ***  engvall  ***
 600  x(1) = 1.0d0
      x(2) = 2.0d0
      x(3) = 0.0d0
      go to 9999
c  *** branin  ***
 700  x(1) = 2.0d0
      x(2) = 0.0d0
      go to 9999
c  ***  beale  ***
 800  x(1) = 1.0d-1
      x(2) = 1.0d-1
      go to 9999
c  *** cragg and levy  ***
 900  x(1) = 1.0d0
      x(2) = 2.0d0
      x(3) = 2.0d0
      x(4) = 2.0d0
      go to 9999
c  ***  box  ***
 1000 x(1) = 0.0d0
      x(2) = 1.0d1
      x(3) = 2.0d1
      go to 9999
c  ***  davidon 1  ***
 1100 do 1101 i = 1,p
 1101    x(i) = 0.0d0
      go to 9999
c  ***  freudenstein and roth  ***
 1200 x(1) = 1.5d1
      x(2) = -2.0d0
      go to 9999
c  ***  watson  ***
 1300 continue
 1400 continue
 1500 continue
 1600 do 1601 i = 1,p
 1601    x(i) = 0.0d0
      go to 9999
c  ***  chebyquad  ***
 1700 pp1inv = 1.0d0/dfloat(p + 1)
      do 1701 i = 1, p
 1701    x(i) = dfloat(i)*pp1inv
      go to 9999
c  *** brown and dennis  ***
 1800 x(1) = 2.5d1
      x(2) = 5.0d0
      x(3) = -5.0d0
      x(4) = -1.0d0
      go to 9999
c  ***  bard  ***
 1900 x(1) = 1.d0
      x(2) = 1.d0
      x(3) = 1.d0
      go to 9999
c  ***  jennrich and sampson  ***
 2000 x(1) = 3.0d-1
      x(2) = 4.0d-1
      go to 9999
c  ***  kowalik and osborne  ***
 2100 x(1) = 2.5d-1
      x(2) = 3.9d-1
      x(3) = 4.15d-1
      x(4) = 3.9d-1
      go to 9999
c  ***  osborne 1  ***
 2200 x(1) = 5.0d-1
      x(2) = 1.5d0
      x(3) = -1.0d0
      x(4) = 1.0d-2
      x(5) = 2.0d-2
      go to 9999
c  ***  osborne 2  ***
 2300 x(1) = 1.3d0
      x(2) = 6.5d-1
      x(3) = 6.5d-1
      x(4) = 7.0d-1
      x(5) = 6.0d-1
      x(6) = 3.0d0
      x(7) = 5.0d0
      x(8) = 7.0d0
      x(9) = 2.0d0
      x(10) = 4.5d0
      x(11) = 5.5d0
      go to 9999
c  ***  madsen  ***
 2400 x(1) = 3.0d0
      x(2) = 1.0d0
      go to 9999
c  ***  meyer  **
 2500 x(1) = 2.0d-2
      x(2) = 4.0d3
      x(3) = 2.5d2
      go to 9999
c  ***  brown  ***
 2600 continue
 2700 continue
 2800 continue
 2900 do 2901 i = 1, p
 2901    x(i) = 5.d-1
      go to 9999
c
c
 9999 return
      end
