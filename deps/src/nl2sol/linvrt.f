      subroutine linvrt(n, lin, l)
c
c  ***  compute  lin = l**-1,  both  n x n  lower triang. stored   ***
c  ***  compactly by rows.  lin and l may share the same storage.  ***
c
c  ***  parameters  ***
c
      integer n
      double precision l(1), lin(1)
c     dimension l(n*(n+1)/2), lin(n*(n+1)/2)
c
c  ***  local variables  ***
c
      integer i, ii, im1, jj, j0, j1, k, k0, np1
      double precision one, t, zero
c/6
c      data one/1.d+0/, zero/0.d+0/
c/7
      parameter (one=1.d+0, zero=0.d+0)
c/
c
c  ***  body  ***
c
      np1 = n + 1
      j0 = n*(np1)/2
      do 30 ii = 1, n
         i = np1 - ii
         lin(j0) = one/l(j0)
         if (i .le. 1) go to 999
         j1 = j0
         im1 = i - 1
         do 20 jj = 1, im1
              t = zero
              j0 = j1
              k0 = j1 - jj
              do 10 k = 1, jj
                   t = t - l(k0)*lin(j0)
                   j0 = j0 - 1
                   k0 = k0 + k - i
 10                continue
              lin(j0) = t/l(k0)
 20           continue
         j0 = j0 - 1
 30      continue
 999  return
c  ***  last card of linvrt follows  ***
      end
