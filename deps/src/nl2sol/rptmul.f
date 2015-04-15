      subroutine rptmul(func, ipivot, j, nn, p, rd, x, y, z)
c
c  ***  func = 1... set  y = rmat * (perm**t) * x.
c  ***  func = 2... set  y = perm * (rmat**t) * rmat * (perm**t) * x.
c  ***  func = 3... set  y = perm * (rmat**t) x.
c
c
c  ***  perm = matrix whose i-th col. is the ipivot(i)-th unit vector.
c  ***  rmat is the upper triangular matrix whose strict upper triangle
c  ***       is stored in  j  and whose diagonal is stored in rd.
c  ***  z is a scratch vector.
c  ***  x and y may share storage.
c
      integer func, nn, p
      integer ipivot(p)
      double precision j(nn,p), rd(p), x(p), y(p), z(p)
c
c  ***  local variables  ***
c
      integer i, im1, k, km1
      double precision zk
c
c  ***  external function  ***
c
      external dotprd
      double precision dotprd
c
c-----------------------------------------------------------------------
c
      if (func .gt. 2) go to 50
c
c  ***  first set  z = (perm**t) * x  ***
c
      do 10 i = 1, p
         k = ipivot(i)
         z(i) = x(k)
 10      continue
c
c  ***  now set  y = rmat * z  ***
c
      y(1) = z(1) * rd(1)
      if (p .le. 1) go to 40
      do 30 k = 2, p
         km1 = k - 1
         zk = z(k)
         do 20 i = 1, km1
 20           y(i) = y(i) + j(i,k)*zk
         y(k) = zk*rd(k)
 30      continue
c
 40   if (func .le. 1) go to 999
      go to 70
c
 50   do 60 i = 1, p
 60      y(i) = x(i)
c
c  ***  set  z = (rmat**t) * y  ***
c
 70   z(1) = y(1) * rd(1)
      if (p .eq. 1) go to 90
      do 80 i = 2, p
         im1 = i - 1
         z(i) = y(i)*rd(i) + dotprd(im1, j(1,i), y)
 80      continue
c
c  ***  now set  y = perm * z  ***
c
 90   do 100 i = 1, p
         k = ipivot(i)
         y(k) = z(i)
 100     continue
c
 999  return
c  ***  last card of rptmul follows  ***
      end
