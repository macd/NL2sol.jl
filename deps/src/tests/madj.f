      subroutine madj(n, p, x, nf, j, uiparm, urparm, ufparm)
      integer n, p, nf, uiparm(1)
      double precision x(p), j(n,p), urparm(1)
      external ufparm
      j(1,1) = 2.0d0*x(1) + x(2)
      j(1,2) = 2.0d0*x(2) + x(1)
      j(2,1) = dcos(x(1))
      j(2,2) = 0.0d0
      j(3,1) = 0.0d0
      j(3,2) = -dsin(x(2))
      return
      end
