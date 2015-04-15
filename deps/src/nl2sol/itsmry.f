      subroutine itsmry(d, iv, p, v, x)
c
c  ***  print nl2sol (version 2.2) iteration summary  ***
c
c  ***  parameter declarations  ***
c
      integer iv(1), p
      double precision d(p), v(1), x(p)
c     dimension iv(*), v(*)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c  ***  local variables  ***
c
      integer cov1, g1, i, ii, iv1, i1, j, m, nf, ng, ol, pu
c/6
c      real model1(6), model2(6)
c/7
      character*4 model1(6), model2(6)
c/
      double precision nreldf, oldf, preldf, reldf, zero
c
c  ***  intrinsic functions  ***
c/+
      integer iabs
c/
c  ***  no external functions or subroutines  ***
c
c  ***  subscripts for iv and v  ***
c
      integer covmat, covprt, covreq, dstnrm, f, fdif, f0, g,
     1        needhd, nfcall, nfcov, ngcov, ngcall, niter, nreduc,
     2        outlev, preduc, prntit, prunit, reldx, size, solprt,
     3        statpr, stppar, sused, x0prt
c
c  ***  iv subscript values  ***
c
c/6
c      data covmat/26/, covprt/14/, g/28/, covreq/15/,
c     1     needhd/39/, nfcall/6/, nfcov/40/, ngcov/41/,
c     2     ngcall/30/, niter/31/, outlev/19/, prntit/48/,
c     3     prunit/21/, solprt/22/, statpr/23/, sused/57/,
c     4     x0prt/24/
c/7
      parameter (covmat=26, covprt=14, g=28, covreq=15,
     1     needhd=39, nfcall=6, nfcov=40, ngcov=41,
     2     ngcall=30, niter=31, outlev=19, prntit=48,
     3     prunit=21, solprt=22, statpr=23, sused=57,
     4     x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
c      data dstnrm/2/, f/10/, f0/13/, fdif/11/, nreduc/6/,
c     1     preduc/7/, reldx/17/, size/47/, stppar/5/
c/7
      parameter (dstnrm=2, f=10, f0=13, fdif=11, nreduc=6,
     1     preduc=7, reldx=17, size=47, stppar=5)
c/
c
c/6
c      data zero/0.d+0/
c/7
      parameter (zero=0.d+0)
c/
c/6
c      data model1(1)/4h    /, model1(2)/4h    /, model1(3)/4h    /,
c     1     model1(4)/4h    /, model1(5)/4h  g /, model1(6)/4h  s /,
c     2     model2(1)/4h g  /, model2(2)/4h s  /, model2(3)/4hg-s /,
c     3     model2(4)/4hs-g /, model2(5)/4h-s-g/, model2(6)/4h-g-s/
c/7
      data model1/'    ','    ','    ','    ','  g ','  s '/,
     1     model2/' g  ',' s  ','g-s ','s-g ','-s-g','-g-s'/
c/
c
c-----------------------------------------------------------------------
c
      pu = iv(prunit)
      if (pu .eq. 0) go to 999
      iv1 = iv(1)
      ol = iv(outlev)
      if (iv1 .lt. 2 .or. iv1 .gt. 15) go to 140
      if (ol .eq. 0) go to 20
      if (iv1 .ge. 12) go to 20
      if (iv1 .ge. 10 .and. iv(prntit) .eq. 0) go to 20
      if (iv1 .gt. 2) go to 10
         iv(prntit) = iv(prntit) + 1
         if (iv(prntit) .lt. iabs(ol)) go to 999
 10   nf = iv(nfcall) - iabs(iv(nfcov))
      iv(prntit) = 0
      reldf = zero
      preldf = zero
      oldf = v(f0)
      if (oldf .le. zero) go to 12
         reldf = v(fdif) / oldf
         preldf = v(preduc) / oldf
 12   if (ol .gt. 0) go to 15
c
c        ***  print short summary line  ***
c
         if (iv(needhd) .eq. 1) write(pu, 1010)
 1010 format(12h0   it    nf,6x,1hf,8x,5hreldf,6x,6hpreldf,5x,5hreldx)
         iv(needhd) = 0
         write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx)
         go to 20
c
c     ***  print long summary line  ***
c
 15   if (iv(needhd) .eq. 1) write(pu,1015)
 1015 format(12h0   it    nf,6x,1hf,8x,5hreldf,6x,6hpreldf,5x,5hreldx,
     1       4x,15hmodel    stppar,6x,4hsize,6x,6hd*step,5x,7hnpreldf)
      iv(needhd) = 0
      m = iv(sused)
      nreldf = zero
      if (oldf .gt. zero) nreldf = v(nreduc) / oldf
      write(pu,1017) iv(niter), nf, v(f), reldf, preldf, v(reldx),
     1               model1(m), model2(m), v(stppar), v(size),
     2               v(dstnrm), nreldf
 1017 format(1x,i5,i6,4d11.3,a3,a4,4d11.3)
c
 20   go to (999,999,30,35,40,45,50,60,70,80,90,150,110,120,130), iv1
c
 30   write(pu,1030)
 1030 format(26h0***** x-convergence *****)
      go to 180
c
 35   write(pu,1035)
 1035 format(42h0***** relative function convergence *****)
      go to 180
c
 40   write(pu,1040)
 1040 format(49h0***** x- and relative function convergence *****)
      go to 180
c
 45   write(pu,1045)
 1045 format(42h0***** absolute function convergence *****)
      go to 180
c
 50   write(pu,1050)
 1050 format(33h0***** singular convergence *****)
      go to 180
c
 60   write(pu,1060)
 1060 format(30h0***** false convergence *****)
      go to 180
c
 70   write(pu,1070)
 1070 format(38h0***** function evaluation limit *****)
      go to 180
c
 80   write(pu,1080)
 1080 format(28h0***** iteration limit *****)
      go to 180
c
 90   write(pu,1090)
 1090 format(18h0***** stopx *****)
      go to 180
c
 110  write(pu,1100)
 1100 format(45h0***** initial sum of squares overflows *****)
c
      go to 150
c
 120  write(pu,1120)
 1120 format(37h0***** bad parameters to assess *****)
      go to 999
c
 130  write(pu,1130)
 1130 format(36h0***** j could not be computed *****)
      if (iv(niter) .gt. 0) go to 190
      go to 150
c
 140  write(pu,1140) iv1
 1140 format(14h0***** iv(1) =,i5,6h *****)
      go to 999
c
c  ***  initial call on itsmry  ***
c
 150  if (iv(x0prt) .ne. 0) write(pu,1150) (i, x(i), d(i), i = 1, p)
 1150 format(23h0    i     initial x(i),7x,4hd(i)//(1x,i5,d17.6,d14.3))
      if (iv1 .ge. 13) go to 999
      iv(needhd) = 0
      iv(prntit) = 0
      if (ol .eq. 0) go to 999
      if (ol .lt. 0) write(pu,1010)
      if (ol .gt. 0) write(pu,1015)
      write(pu,1160) v(f)
 1160 format(12h0    0     1,d11.3,11x,d11.3)
      go to 999
c
c  ***  print various information requested on solution  ***
c
 180  iv(needhd) = 1
      if (iv(statpr) .eq. 0) go to 190
         oldf = v(f0)
         preldf = zero
         nreldf = zero
         if (oldf .le. zero) go to 185
              preldf = v(preduc) / oldf
              nreldf = v(nreduc) / oldf
 185     nf = iv(nfcall) - iv(nfcov)
         ng = iv(ngcall) - iv(ngcov)
         write(pu,1180) v(f), v(reldx), nf, ng, preldf, nreldf
 1180 format(9h0function,d17.6,8h   reldx,d20.6/12h func. evals,
     1   i8,9x,11hgrad. evals,i8/7h preldf,d19.6,3x,7hnpreldf,d18.6)
c
         if (iv(nfcov) .gt. 0) write(pu,1185) iv(nfcov)
 1185    format(1h0,i4,34h extra func. evals for covariance.)
         if (iv(ngcov) .gt. 0) write(pu,1186) iv(ngcov)
 1186    format(1x,i4,34h extra grad. evals for covariance.)
c
 190  if (iv(solprt) .eq. 0) go to 210
         iv(needhd) = 1
         g1 = iv(g)
         write(pu,1190)
 1190 format(22h0    i      final x(i),8x,4hd(i),10x,4hg(i)/)
         do 200 i = 1, p
              write(pu,1200) i, x(i), d(i), v(g1)
              g1 = g1 + 1
 200          continue
 1200    format(1x,i5,d17.6,2d14.3)
c
 210  if (iv(covprt) .eq. 0) go to 999
      cov1 = iv(covmat)
      iv(needhd) = 1
      if (cov1) 220, 230, 240
 220  if (-1 .eq. cov1) write(pu,1220)
 1220 format(43h0++++++ indefinite covariance matrix ++++++)
      if (-2 .eq. cov1) write(pu,1225)
 1225 format(52h0++++++ oversize steps in computing covariance +++++)
      go to 999
c
 230  write(pu,1230)
 1230 format(45h0++++++ covariance matrix not computed ++++++)
      go to 999
c
 240  i = iabs(iv(covreq))
      if (i .le. 1) write(pu,1241)
 1241 format(48h0covariance = scale * h**-1 * (j**t * j) * h**-1/)
      if (i .eq. 2) write(pu,1242)
 1242 format(27h0covariance = scale * h**-1/)
      if (i .ge. 3) write(pu,1243)
 1243 format(36h0covariance = scale * (j**t * j)**-1/)
      ii = cov1 - 1
      if (ol .le. 0) go to 260
      do 250 i = 1, p
         i1 = ii + 1
         ii = ii + i
         write(pu,1250) i, (v(j), j = i1, ii)
 250     continue
 1250 format(4h row,i3,2x,9d12.4/(9x,9d12.4))
      go to 999
c
 260  do 270 i = 1, p
         i1 = ii + 1
         ii = ii + i
         write(pu,1270) i, (v(j), j = i1, ii)
 270     continue
 1270 format(4h row,i3,2x,5d12.4/(9x,5d12.4))
c
 999  return
c  ***  last card of itsmry follows  ***
      end
