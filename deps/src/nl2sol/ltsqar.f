      subroutine ltsqar(n, a, l)
c
c  ***  set a to lower triangle of (l**t) * l  ***
c
c  ***  l = n x n lower triang. matrix stored rowwise.  ***
c  ***  a is also stored rowwise and may share storage with l.  ***
c
      integer n
      double precision a(1), l(1)
c     dimension a(n*(n+1)/2), l(n*(n+1)/2)
c
      integer i, ii, iim1, i1, j, k, m
      double precision lii, lj
c
      ii = 0
      do 50 i = 1, n
         i1 = ii + 1
         ii = ii + i
         m = 1
         if (i .eq. 1) go to 30
         iim1 = ii - 1
         do 20 j = i1, iim1
              lj = l(j)
              do 10 k = i1, j
                   a(m) = a(m) + lj*l(k)
                   m = m + 1
 10                continue
 20           continue
 30      lii = l(ii)
         do 40 j = i1, ii
 40           a(j) = lii * l(j)
 50      continue
c
 999  return
c  ***  last card of ltsqar follows  ***
      end
