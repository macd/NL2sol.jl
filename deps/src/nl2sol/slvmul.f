      subroutine slvmul(p, y, s, x)
c
c  ***  set  y = s * x,  s = p x p symmetric matrix.  ***
c  ***  lower triangle of  s  stored rowwise.         ***
c
c  ***  parameter declarations  ***
c
      integer p
      double precision s(1), x(p), y(p)
c     dimension s(p*(p+1)/2)
c
c  ***  local variables  ***
c
      integer i, im1, j, k
      double precision xi
c
c  ***  no intrinsic functions  ***
c
c  ***  external function  ***
c
      external dotprd
      double precision dotprd
c
c-----------------------------------------------------------------------
c
      j = 1
      do 10 i = 1, p
         y(i) = dotprd(i, s(j), x)
         j = j + i
 10      continue
c
      if (p .le. 1) go to 999
      j = 1
      do 40 i = 2, p
         xi = x(i)
         im1 = i - 1
         j = j + 1
         do 30 k = 1, im1
              y(k) = y(k) + s(j)*xi
              j = j + 1
 30           continue
 40      continue
c
 999  return
c  ***  last card of slvmul follows  ***
      end
