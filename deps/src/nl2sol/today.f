      subroutine today(datime)
c
c  ***  supply sumsol version  ***
c
c/6
c      real datime(4), dt1, dt2, dt3, dt4
c      data dt1,dt2,dt3,dt4/4hnl2s,4hol  ,4hver.,4h2.2 /
c/7
      character*4 datime(4), dt1, dt2, dt3, dt4
      data dt1,dt2,dt3,dt4/'nl2s','ol  ','ver.','2.2 '/
c/
c
      datime(1) = dt1
      datime(2) = dt2
      datime(3) = dt3
      datime(4) = dt4
 999  return
c  ***  last line of datime follows  ***
      end
