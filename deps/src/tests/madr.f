      subroutine madr(n, p, x, nf, r, uiparm, urparm, ufparm)
      integer n, p, nf, uiparm(1)
      double precision x(p), r(n), urparm(1)
      external ufparm
      r(1) = x(1)**2 + x(2)**2 + x(1)*x(2)
      r(2) = dsin(x(1))
      r(3) = dcos(x(2))
      return
      end
