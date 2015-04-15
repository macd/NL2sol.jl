      subroutine nltest (n, p, nex, title1, title2, rstart)
c
c  ***  call nl2sol, save and print statistics  ***
c
c
      integer n, p, nex
      logical rstart
c/6
c      real title1, title2
c/7
      character*4 title1, title2
c/
c
      common /testcm/ v, rs, jac, nout, nprob, xscal1, xscal2, is, iv
      common /testch/ name, irc
      integer is(6,82), iv(100), jac, nout, nprob, xscal1, xscal2
      real rs(5,82)
c/6
c      integer irc(50)
c      real name(2,50)
c/7
      character name(2,82)*4, irc(82)*1
c/
      double precision v(3000)
c
      logical rstrt
      integer i, irun, pu, uip(1)
c/6
c      integer alg(2), jtyp(2), rc(10)
c      real datime(4)
c/7
      character*4 datime(4)
      character*2 alg(2)
      character*1 jtyp(2), rc(10)
c/
      double precision one, t, urparm(1), x(20), x0scal, zero
c
c     ***  external functions and subroutines  ***
c
      external nl2sno, nl2sol, testr, testj, today, xinit
c
c  ***  iv and v subscripts  ***
c
      integer f, f0, nfcall, nfcov, ngcall, niter, nreduc, preduc,
     1        prunit, reldx
c
c/6
c      data f/10/, f0/13/, nfcall/6/, nfcov/40/, ngcall/30/,
c     1     ngcov/41/, niter/31/, nreduc/6/, preduc/7/,
c     2     prunit/21/, reldx/17/
c/7
      parameter (f=10, f0=13, nfcall=6, nfcov=40, ngcall=30,
     1     ngcov=41, niter=31, nreduc=6, preduc=7,
     2     prunit=21, reldx=17)
c/
c/6
c      data one/1.d+0/, zero/0.d+0/
c/7
      parameter (one=1.d+0, zero=0.d+0)
c/
c/6
c      data alg(1),alg(2)/2hol,2hno/, jtyp(1),jtyp(2)/1h ,1h*/
c      data rc(1)/1h./, rc(2)/1h+/, rc(3)/1hx/, rc(4)/1hr/, rc(5)/1hb/,
c     1     rc(6)/1ha/, rc(7)/1hs/, rc(8)/1hf/, rc(9)/1he/, rc(10)/1hi/
c/7
      data alg(1),alg(2)/'ol','no'/, jtyp(1),jtyp(2)/' ','*'/
      data rc(1)/'.'/, rc(2)/'+'/, rc(3)/'x'/, rc(4)/'r'/, rc(5)/'b'/,
     1     rc(6)/'a'/, rc(7)/'s'/, rc(8)/'f'/, rc(9)/'e'/, rc(10)/'i'/
c/
c
c-----------------------------------------------------------------------
c
      uip(1) = nex
      rstrt = rstart
 160  format(a15)
      if (rstrt) go to 20
         pu = iv(prunit)
         call today(datime)
         if (pu .ne. 0) write(pu,10) alg(jac), title1, title2, datime
 10      format (1h1//11h ***** nl2s,a2,12h on problem ,2a4,6h *****,6x,
     1           2a4,2x,2a4)
c
 20   do 100 irun = xscal1, xscal2
         if (rstrt) go to 40
         iv(1) = 12
         x0scal = 1.0d1 ** (irun-1)
c
c        ***  initialize the solution vector x  ***
         call xinit(p, x, nex)
 180     format(a21, i5)
         do 30 i = 1, p
 30           x(i) = x0scal * x(i)
c
 40      if (jac .eq. 1)
     1             call nl2sol(n,p,x,testr,testj,iv,v,uip,urparm,testr)
         if (jac .eq. 2)
     1             call nl2sno(n,p,x,testr,iv,v,uip,urparm,testr)
         if (.not. rstrt .and. nprob .lt. 82) nprob = nprob + 1
         name(1,nprob) = title1
         name(2,nprob) = title2
         is(1,nprob) = n
         is(2,nprob) = p
         is(3,nprob) = iv(niter)
         is(4,nprob) = iv(nfcall) - iv(nfcov)
         is(5,nprob) = iv(ngcall) - iv(ngcov)
         i = iv(1)
         irc(nprob) = rc(i)
         is(6,nprob) = jac
         rs(1,nprob) = x0scal
         rs(2,nprob) = v(f)
         t = one
         if (v(f0) .gt. zero) t = v(preduc) / v(f0)
         rs(3,nprob) = t
         t = one
         if (v(f0) .gt. zero) t = v(nreduc) / v(f0)
         rs(4,nprob) = t
         rs(5,nprob) = v(reldx)
         rstrt = .false.
         if (nout .eq. 0) go to 100
         if (nprob .eq. 1) write(nout,50) datime
 50      format(1h1,11x,2a4,2x,2a4,10x,24hnl2sol test summary.....,10x,
     1          32h(* = finite-difference jacobian)/
     2          48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x,
     3          39hfinal f     preldf     nreldf     reldx/)
         write(nout,60) jtyp(jac), title1, title2,
     1                (is(i,nprob),i=1,5),irc(nprob),(rs(i,nprob),i=1,5)
 60      format(1x,a1,2a4,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)
 100     continue
c
 999  return
c  ***  last card of nltest follows  ***
      end
