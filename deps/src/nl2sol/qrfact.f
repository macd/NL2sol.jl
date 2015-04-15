      subroutine qrfact(nm,m,n,qr,alpha,ipivot,ierr,nopivk,sum)
c
c  ***  compute the qr decomposition of the matrix stored in qr  ***
c
c     *****parameters.
      integer nm,m,n,ipivot(n),ierr,nopivk
      double precision  qr(nm,n),alpha(n),sum(n)
c     *****local variables.
      integer i,j,jbar,k,k1,minum,mk1
      double precision  alphak,beta,qrkk,qrkmax,sigma,temp,ufeta,rktol,
     1        rktol1,sumj
c     *****functions.
c/+
      integer min0
      double precision  dabs,dsqrt
c/
      external dotprd, rmdcon, vaxpy, vscopy, v2norm
      double precision dotprd, rmdcon, v2norm
c dotprd... returns inner product of two vectors.
c rmdcon... returns machine-dependent constants.
c vaxpy... computes scalar times one vector plus another.
c vscopy... sets all elements of a vector to a scalar.
c v2norm... returns the 2-norm of a vector.
c
c     *****constants.
      double precision one, p01, p99, zero
c/6
c      data one/1.0d+0/, p01/0.01d+0/, p99/0.99d+0/, zero/0.0d+0/
c/7
      parameter (one=1.0d+0, p01=0.01d+0, p99=0.99d+0, zero=0.0d+0)
      save rktol, ufeta
c/
c
c
c     ..................................................................
c     ..................................................................
c
c
c     *****purpose.
c
c     this subroutine does a qr-decomposition on the m x n matrix qr,
c        with an optionally modified column pivoting, and returns the
c        upper triangular r-matrix, as well as the orthogonal vectors
c        used in the transformations.
c
c     *****parameter description.
c     on input.
c
c        nm must be set to the row dimension of the two dimensional
c             array parameters as declared in the calling program
c             dimension statement.
c
c        m must be set to the number of rows in the matrix.
c
c        n must be set to the number of columns in the matrix.
c
c        qr contains the real rectangular matrix to be decomposed.
c
c     nopivk is used to control pivotting.  columns 1 through
c        nopivk will remain fixed in position.
c
c        sum is used for temporary storage for the subroutine.
c
c     on output.
c
c        qr contains the non-diagonal elements of the r-matrix
c             in the strict upper triangle. the vectors u, which
c             define the householder transformations   i - u*u-transp,
c             are in the columns of the lower triangle. these vectors u
c             are scaled so that the square of their 2-norm is 2.0.
c
c        alpha contains the diagonal elements of the r-matrix.
c
c        ipivot reflects the column pivoting performed on the input
c             matrix to accomplish the decomposition. the j-th
c             element of ipivot gives the column of the original
c             matrix which was pivoted into column j during the
c             decomposition.
c
c        ierr is set to.
c             0 for normal return,
c             k if no non-zero pivot could be found for the k-th
c                  transformation, or
c             -k for an error exit on the k-th thansformation.
c             if an error exit was taken, the first (k - 1)
c             transformations are correct.
c
c
c     *****applications and usage restrictions.
c     this may be used when solving linear least-squares problems --
c     see subroutine qr1 of rosepack.  it is called for this purpose
c     by llsqst in the nl2sol (nonlinear least-squares) package.
c
c     *****algorithm notes.
c     this version of qrfact tries to eliminate the occurrence of
c     underflows during the accumulation of inner products.  rktol1
c     is chosen below so as to insure that discarded terms have no
c     effect on the computed two-norms.
c
c     adapted from the algol routine solve (1).
c
c     *****references.
c     (1)     businger,p. and golub,g.h., linear least squares
c     solutions by housholder transformations, in wilkinson,j.h.
c     and reinsch,c.(eds.), handbook for automatic computation,
c     volume ii. linear algebra, springer-verlag, 111-118 (1971).
c     prepublished in numer.math. 7, 269-276 (1965).
c
c     *****history.
c     this amounts to the subroutine qr1 of rosepack with rktol1 used
c     in place of rktol below, with v2norm used to initialize (and
c     sometimes update) the sum array, and with calls on dotprd and
c     vaxpy in place of some loops.
c
c     *****general.
c
c     development of this program supported in part by
c     national science foundation grant gj-1154x3 and
c     national science foundation grant dcr75-08802
c     to national bureau of economic research, inc.
c
c
c
c     ..................................................................
c     ..................................................................
c
c
c     ..........  ufeta is the smallest positive floating point number
c        s.t. ufeta and -ufeta can both be represented.
c
c     ..........  rktol is the square root of the relative precision
c        of floating point arithmetic (machep).
      data rktol/0.d+0/, ufeta/0.d+0/
c     *****body of program.
      if (ufeta .gt. zero) go to 10
         ufeta = rmdcon(1)
         rktol = rmdcon(4)
   10 ierr = 0
      rktol1 = p01 * rktol
c
      do 20 j=1,n
         sum(j) = v2norm(m, qr(1,j))
         ipivot(j) = j
   20 continue
c
      minum = min0(m,n)
c
      do 120 k=1,minum
         mk1 = m - k + 1
c        ..........k-th householder transformation..........
         sigma = zero
         jbar = 0
c        ..........find largest column sum..........
      if (k .le. nopivk) go to 50
         do 30 j=k,n
              if (sigma .ge. sum(j))  go to 30
              sigma = sum(j)
              jbar = j
   30    continue
c
         if (jbar .eq. 0)  go to 220
         if (jbar .eq. k)  go to 50
c        ..........column interchange..........
         i = ipivot(k)
         ipivot(k) = ipivot(jbar)
         ipivot(jbar) = i
         sum(jbar) = sum(k)
         sum(k) = sigma
c
         do 40 i=1,m
              sigma = qr(i,k)
              qr(i,k) = qr(i,jbar)
              qr(i,jbar) = sigma
   40    continue
c        ..........end of column interchange..........
   50    continue
c        ..........  second inner product  ..........
         qrkmax = zero
c
         do 60 i=k,m
              if (dabs( qr(i,k) ) .gt. qrkmax)  qrkmax = dabs( qr(i,k) )
   60    continue
c
         if (qrkmax .lt. ufeta)  go to 210
         alphak = v2norm(mk1, qr(k,k)) / qrkmax
         sigma = alphak**2
c
c        ..........  end second inner product  ..........
         qrkk = qr(k,k)
         if (qrkk .ge. zero)  alphak = -alphak
         alpha(k) = alphak * qrkmax
         beta = qrkmax * dsqrt(sigma - (qrkk*alphak/qrkmax) )
         qr(k,k) = qrkk - alpha(k)
         do 65 i=k,m
   65         qr(i,k) =  qr(i,k) / beta
         k1 = k + 1
         if (k1 .gt. n) go to 120
c
         do 110 j = k1, n
              temp = -dotprd(mk1, qr(k,k), qr(k,j))
c
c             ***  set qr(i,j) = qr(i,j) + temp*qr(i,k), i = k,...,m.
c
              call vaxpy(mk1, qr(k,j), temp, qr(k,k), qr(k,j))
c
              if (k1 .gt. m) go to 110
              sumj = sum(j)
              if (sumj .lt. ufeta) go to 110
              temp = dabs(qr(k,j)/sumj)
              if (temp .lt. rktol1) go to 110
              if (temp .ge. p99) go to 90
                   sum(j) = sumj * dsqrt(one - temp**2)
                   go to 110
   90         sum(j) = v2norm(m-k, qr(k1,j))
  110    continue
c        ..........end of k-th householder transformation..........
  120 continue
c
      go to 999
c     ..........error exit on k-th transformation..........
  210 ierr = -k
      go to 230
c     ..........no non-zero acceptable pivot found..........
  220 ierr = k
  230 do 240 i = k, n
         alpha(i) = zero
         if (i .gt. k) call vscopy(i-k, qr(k,i), zero)
 240     continue
c     ..........return to caller..........
  999 return
c     ..........last card of qrfact..........
      end
