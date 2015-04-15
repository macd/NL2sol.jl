      subroutine dfault(iv, v)
c
c  ***  supply nl2sol (version 2.2) default values to iv and v  ***
c
      integer iv(25)
      double precision v(45)
c/+
      double precision dmax1
c/
      external imdcon, rmdcon
      integer imdcon
      double precision rmdcon
c
      double precision machep, mepcrt, one, sqteps, three
c
c  ***  subscripts for iv and v  ***
c
      integer afctol, cosmin, covprt, covreq, decfac, delta0, dfac,
     1        dinit, dltfdc, dltfdj, dtype, d0init, epslon, fuzz,
     2        incfac, inits, jtinit, lmax0, mxfcal, mxiter, outlev,
     3        parprt, phmnfc, phmxfc, prunit, rdfcmn, rdfcmx,
     4        rfctol, rlimit, solprt, statpr, tuner1, tuner2, tuner3,
     5        tuner4, tuner5, xctol, xftol, x0prt
c
c/6
c      data one/1.d+0/, three/3.d+0/
c/7
      parameter (one=1.d+0, three=3.d+0)
c/
c
c  ***  iv subscript values  ***
c
c/6
c      data covprt/14/, covreq/15/, dtype/16/, inits/25/,
c     1     mxfcal/17/, mxiter/18/, outlev/19/,
c     2     parprt/20/, prunit/21/, solprt/22/,
c     3     statpr/23/, x0prt/24/
c/7
      parameter (covprt=14, covreq=15, dtype=16, inits=25,
     1     mxfcal=17, mxiter=18, outlev=19,
     2     parprt=20, prunit=21, solprt=22,
     3     statpr=23, x0prt=24)
c/
c
c  ***  v subscript values  ***
c
c/6
c      data afctol/31/, cosmin/43/, decfac/22/,
c     1     delta0/44/, dfac/41/, dinit/38/, dltfdc/40/,
c     2     dltfdj/36/, d0init/37/, epslon/19/, fuzz/45/,
c     3     incfac/23/, jtinit/39/, lmax0/35/, phmnfc/20/,
c     4     phmxfc/21/, rdfcmn/24/, rdfcmx/25/,
c     5     rfctol/32/, rlimit/42/, tuner1/26/,
c     6     tuner2/27/, tuner3/28/, tuner4/29/,
c     7     tuner5/30/, xctol/33/, xftol/34/
c/7
      parameter (afctol=31, cosmin=43, decfac=22,
     1     delta0=44, dfac=41, dinit=38, dltfdc=40,
     2     dltfdj=36, d0init=37, epslon=19, fuzz=45,
     3     incfac=23, jtinit=39, lmax0=35, phmnfc=20,
     4     phmxfc=21, rdfcmn=24, rdfcmx=25,
     5     rfctol=32, rlimit=42, tuner1=26,
     6     tuner2=27, tuner3=28, tuner4=29,
     7     tuner5=30, xctol=33, xftol=34)
c/
c
c-----------------------------------------------------------------------
c
      iv(1) = 12
      iv(covprt) = 1
      iv(covreq) = 1
      iv(dtype) = 1
      iv(inits) = 0
      iv(mxfcal) = 200
      iv(mxiter) = 150
      iv(outlev) = 1
      iv(parprt) = 1
      iv(prunit) = imdcon(1)
      iv(solprt) = 1
      iv(statpr) = 1
      iv(x0prt) = 1
c
      machep = rmdcon(3)
      v(afctol) = 1.d-20
      if (machep .gt. 1.d-10) v(afctol) = machep**2
      v(cosmin) = dmax1(1.d-6, 1.d+2 * machep)
      v(decfac) = 0.5d+0
      sqteps = rmdcon(4)
      v(delta0) = sqteps
      v(dfac) = 0.6d+0
      v(dinit) = 0.d+0
      mepcrt = machep ** (one/three)
      v(dltfdc) = mepcrt
      v(dltfdj) = sqteps
      v(d0init) = 1.d+0
      v(epslon) = 0.1d+0
      v(fuzz) = 1.5d+0
      v(incfac) = 2.d+0
      v(jtinit) = 1.d-6
      v(lmax0) = 100.d+0
      v(phmnfc) = -0.1d+0
      v(phmxfc) = 0.1d+0
      v(rdfcmn) = 0.1d+0
      v(rdfcmx) = 4.d+0
      v(rfctol) = dmax1(1.d-10, mepcrt**2)
      v(rlimit) = rmdcon(5)
      v(tuner1) = 0.1d+0
      v(tuner2) = 1.d-4
      v(tuner3) = 0.75d+0
      v(tuner4) = 0.5d+0
      v(tuner5) = 0.75d+0
      v(xctol) = sqteps
      v(xftol) = 1.d+2 * machep
c
 999  return
c  ***  last card of dfault follows  ***
      end
