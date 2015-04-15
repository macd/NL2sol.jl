      subroutine qapply(nn, n, p, j, r, ierr)
c     *****parameters.
      integer nn, n, p, ierr
      double precision j(nn,p), r(n)
c
c     ..................................................................
c     ..................................................................
c
c     *****purpose.
c     this subroutine applies to r the orthogonal transformations
c     stored in j by qrfact
c
c     *****parameter description.
c     on input.
c
c        nn is the row dimension of the matrix j as declared in
c             the calling program dimension statement
c
c        n is the number of rows of j and the size of the vector r
c
c        p is the number of columns of j and the size of sigma
c
c        j contains on and below its diagonal the column vectors
c             u which determine the householder transformations
c             ident - u*u.transpose
c
c        r is the right hand side vector to which the orthogonal
c             transformations will be applied
c
c        ierr if non-zero indicates that not all the transformations
c             were successfully determined and only the first
c             abs(ierr) - 1 transformations will be used
c
c     on output.
c
c        r has been overwritten by its transformed image
c
c     *****application and usage restrictions.
c     none
c
c     *****algorithm notes.
c     the vectors u which determine the householder transformations
c     are normalized so that their 2-norm squared is 2.  the use of
c     these transformations here is in the spirit of (1).
c
c     *****subroutines and functions called.
c
c     dotprd - function, returns the inner product of vectors
c
c     *****references.
c     (1) businger, p. a., and golub, g. h. (1965), linear least squares
c        solutions by householder transformations, numer. math. 7,
c        pp. 269-276.
c
c     *****history.
c     designed by david m. gay, coded by stephen c. peters (winter 1977)
c
c     *****general.
c
c     this subroutine was written in connection with research
c     supported by the national science foundation under grants
c     mcs-7600324, dcr75-10143, 76-14311dss, and mcs76-11989.
c
c     ..................................................................
c     ..................................................................
c
c     *****local variables.
      integer i, k, l, nl1
      double precision t
c     *****intrinsic functions.
c/+
      integer iabs
c/
c     *****functions.
      external dotprd
      double precision dotprd
c
      k = p
      if (ierr .ne. 0) k = iabs(ierr) - 1
      if ( k .eq. 0) go to 999
c
      do 20 l = 1, k
         nl1 = n - l + 1
         t = -dotprd(nl1, j(l,l), r(l))
c
         do 10 i = l, n
 10           r(i) = r(i) + t*j(i,l)
 20   continue
 999  return
c     .... last card of qapply .........................................
      end
