      program nlmain
c///////////////////////////////////////////////////////////////////////
c  ***  run nl2sol on various test problems, print summary statistics.
c
c     *****common storage with nltest.
c
      common /testcm/ v, rs, jac, nout, nprob, xscal1, xscal2, is, iv
      common /testch/ name, irc
      integer is(6,82), iv(100), jac, nout, nprob, xscal1, xscal2
      real rs(5,82)
c/6
c      real name(2,50)
c      integer irc(50)
c/7
      character name(2,82)*4, irc(82)*1
c/
c      double precision v(1736)
      double precision v(3000)
c
c
c     ..................................................................
c
c     *****purpose.
c        this main program calls nltest to run nl2sol, the nonlinear
c     least-squares solver of ref. 1, on various test problems.
c
c
c     *****application and usage restrictions.
c     this main driver is intended to check whether the nl2sol
c     (nonlinear least-squares) package was successfully
c     transported to a new machine.
c
c     *****algorithm notes.
c     the test problems used are from references (2), (3), and (4).
c     some additional test problems were suggested by jorge more (pri-
c     vate communication).  calls passing these problems to nltest have
c     been commented out (since there are enough other problems), but
c     not removed, since they may be of interest to other researchers.
c
c     *****functions and subroutines called.
c
c        dfault - establishes the default parameter settings for
c                 iv and v.
c
c        imdcon - imdcon(2) returns i/o unit number on which nltest
c                  writes a summary of each test run.
c
c        ivvset - supplies nondefault values for iv and v.
c
c        nltest - calls nl2sol, the nonlinear least-squares
c                  problem solver.
c
c        today  - supplies date and time (or current version of nl2sol).
c
c     *****references.
c
c     (1). dennis, j.e.. gay, d.m.. and welsch, r.e. (1980),
c          an adaptive nonlinear least-squares algorithm,
c          submitted to acm trans. math. software.
c          under revision.
c
c     (2). gill, p.e.. and murray, w. (1976),algorithms for the
c          solution of the non-linear least-squares problem,
c          npl report nac71,(national physical laboratory,
c          division of numerical analysis and computing,
c          teddington,middlesex,england).
c
c     (3) meyer, r.r. (1970), theoretical and computational aspects
c        of nonlinear regression, pp. 465-486 of nonlinear programming,
c        edited by j.b. rosen, o.l.mangasarian, and k. ritter,
c        academic press, new york.
c
c     (4) brown, k.m. (1969), a quadratically convergent newton-
c        like method based upon gaussian elimination,
c        siam j. numer. anal. 6, pp. 560-569.
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
c     *****intrinsic functions.
c/+
      integer mod
      real amax1
c/
c     *****external functions and subroutines.
      external dfault, imdcon, ivvset, nltest, today
      integer imdcon
c
c     *****local variables.
      logical rstart
      integer i, j, k, mxfcsv, mxitsv, pu
c/6
c      integer jtyp(2)
c      real datime(4)
c/7
      character datime(4)*4, jtyp(2)*1
c/
c
c/6
c      data rstart/.false./, jtyp(1),jtyp(2)/1h ,1h*/
c/7
      data rstart/.false./, jtyp(1),jtyp(2)/' ','*'/
c/
c
c-----------------------------------------------------------------------
c
c  ***  establish default parameter settings  ***
      call dfault (iv, v)
      nout = imdcon(2)
c
c  ***  non-default parameter settings  ***
c
c      call ivvset(iv, v)
      pu = iv(21)
c
      jac = 1
      nprob = 0
      xscal1 = 1
      xscal2 = 3

c
c/6
c      call nltest(2,2,1,4hrosn,4hbrok,rstart)
c      call nltest(3,3,2,4hheli,4hx   ,rstart)
c      call nltest(4,4,3,4hsing,4hular,rstart)
c      call nltest(7,4,4,4hwood,4hs   ,rstart)
c      xscal2 = 1
c      call nltest(3,3,5,4hzang,4hwill,rstart)
c      xscal2 = 3
c      call nltest(5,3,6,4hengv,4hall ,rstart)
c      call nltest(2,2,7,4hbran,4hin  ,rstart)
c      xscal2 = 2
c      call nltest(3,2,8,4hbeal,4he   ,rstart)
c      call nltest(5,4,9,4hcrag,4hg   ,rstart)
c      xscal2 = 2
c      call nltest(10,3,10,4hbox ,4h    ,rstart)
c      mxfcsv = iv(17)
c      mxitsv = iv(18)
c      iv(17) = 20
c      iv(18) = 15
c      xscal2 = 1
c      call nltest(15,15,11,4hdavi,4hdon1,rstart)
c      iv(17) = mxfcsv
c      iv(18) = mxitsv
c      xscal2 = 3
c      call nltest(2,2,12,4hfrds,4htein,rstart)
c      xscal2 = 1
c      call nltest(31,6,13,4hwats,4hon6 ,rstart)
c      call nltest(31,9,14,4hwats,4hon9 ,rstart)
c      call nltest(31,12,15,4hwats,4hon12,rstart)
c      mxfcsv = iv(17)
c      iv(17) = 20
c      mxitsv = iv(18)
c      iv(18) = 15
c      call nltest(31,20,16,4hwats,4hon20,rstart)
c      iv(17) = mxfcsv
c      iv(18) = mxitsv
c      xscal2 = 2
c      call nltest(8,8,17,4hcheb,4hqd8 ,rstart)
c      xscal2 = 3
c      call nltest(20,4,18,4hbrow,4hn   ,rstart)
c      call nltest(15,3,19,4hbard,4h    ,rstart)
c      xscal2 = 1
c      call nltest(10,2,20,4hjenn,4hrich,rstart)
c      xscal2 = 3
c      call nltest(11,4,21,4hkowa,4hlik ,rstart)
c      xscal2 = 1
c      call nltest(33,5,22,4hosbo,4hrne1,rstart)
c      xscal2 = 2
c      call nltest(65,11,23,4hosbo,4hrne2,rstart)
c      xscal2 = 3
c      call nltest(3,2,24,4hmads,4hen  ,rstart)
c      xscal2 = 1
c      iv(17) = 400
c      iv(18) = 300
c      call nltest(16,3,25,4hmeye,4hr   ,rstart)
c/7
      call nltest(2,2,1,'rosn','brok',rstart)
      call nltest(3,3,2,'heli','x   ',rstart)
      call nltest(4,4,3,'sing','ular',rstart)
      call nltest(7,4,4,'wood','s   ',rstart)
      xscal2 = 1
      call nltest(3,3,5,'zang','will',rstart)
      xscal2 = 3
      call nltest(5,3,6,'engv','all ',rstart)
      call nltest(2,2,7,'bran','in  ',rstart)
      xscal2 = 2
      call nltest(3,2,8,'beal','e   ',rstart)
      call nltest(5,4,9,'crag','g   ',rstart)
      xscal2 = 2
      call nltest(10,3,10,'box ','    ',rstart)
      mxfcsv = iv(17)
      mxitsv = iv(18)
      iv(17) = 20
      iv(18) = 15
      xscal2 = 1
      call nltest(15,15,11,'davi','don1',rstart)
      iv(17) = mxfcsv
      iv(18) = mxitsv
      xscal2 = 3
      call nltest(2,2,12,'frds','tein',rstart)
      xscal2 = 1
      call nltest(31,6,13,'wats','on6 ',rstart)
      call nltest(31,9,14,'wats','on9 ',rstart)
      call nltest(31,12,15,'wats','on12',rstart)
      mxfcsv = iv(17)
      iv(17) = 20
      mxitsv = iv(18)
      iv(18) = 15
      call nltest(31,20,16,'wats','on20',rstart)
      iv(17) = mxfcsv
      iv(18) = mxitsv
      xscal2 = 2
      call nltest(8,8,17,'cheb','qd8 ',rstart)
      xscal2 = 3
      call nltest(20,4,18,'brow','n   ',rstart)
      call nltest(15,3,19,'bard','    ',rstart)
      xscal2 = 1
      call nltest(10,2,20,'jenn','rich',rstart)
      xscal2 = 3
      call nltest(11,4,21,'kowa','lik ',rstart)
      xscal2 = 1
      call nltest(33,5,22,'osbo','rne1',rstart)
      xscal2 = 2
      call nltest(65,11,23,'osbo','rne2',rstart)
      xscal2 = 3
      call nltest(3,2,24,'mads','en  ',rstart)
      xscal2 = 1
      iv(17) = 400
      iv(18) = 300
      call nltest(16,3,25,'meye','r   ',rstart)
c/
c  ***  brown5  ***
c     call nltest(5,5,26,4hbrow,4hn5  ,rstart)
c  ***  brown10  ***
c     call nltest(10,10,27,4hbrow,4hn10 ,rstart)
c  ***  brown30  ***
c     call nltest(30,30,28,4hbrow,4hn30 ,rstart)
c  ***  brown40  ***
c     call nltest(40,40,29,4hbrow,4hn40 ,rstart)
c  ***  bard+10 ***
c     call nltest(15,3,30,4hbard,4h+10 ,rstart)
c  ***  kowalik and osborne + 10  ***
c     call nltest(11,4,31,4hkowa,4hl+10,rstart)
c  ***  meyer + 10  ***
c     call nltest(16,3,32,4hmeye,4hr+10,rstart)
c  ***  watson6 + 10  ***
c     call nltest(31,6,33,4hwat6,4h+10 ,rstart)
c  ***  watson9 + 10  ***
c     call nltest(31,9,34,4hwat9,4h+10 ,rstart)
c  ***  watson12 + 10  ***
c     call nltest(31,12,35,4hwat1,4h2+10,rstart)
c  ***  watson20 + 10  ***
c     call nltest(31,20,36,4hwat2,4h0+10,rstart)



c  ***  brown5  ***
      call nltest(5, 5, 26,'brow','n5  ',rstart)
c  ***  brown10  ***
      call nltest(10,10,27,'brow','n10 ',rstart)
c  ***  brown30  ***
      call nltest(30,30,28,'brow','n30 ',rstart)
c  ***  brown40  ***
c      call nltest(40,40,29,'brow','n40 ',rstart)
c  ***  bard+10 ***
c     call nltest(15,3,30,'bard','+10 ',rstart)
c  ***  kowalik and osborne + 10  ***
c       call nltest(11,4,31,'kowa','l+10',rstart)
c  ***  meyer + 10  ***
c       call nltest(16,3,32,'meye','r+10',rstart)
c  ***  watson6 + 10  ***
c       call nltest(31,6,33,'wat6','+10 ',rstart)
c  ***  watson9 + 10  ***
c       call nltest(31,9,34,'wat9','+10 ',rstart)
c  ***  watson12 + 10  ***
c       call nltest(31,12,35,'wat1','2+10',rstart)
c  ***  watson20 + 10  ***
c       call nltest(31,20,36,'wat2','0+10',rstart)
c



c
c  ***  repeat two tests using finite-difference jacobian  ***
c
      jac = 2
      xscal2 = 1
c
      iv(17) = 50
      iv(18) = 40
c/6
c      call nltest(2,2,1,4hrosn,4hbrok,rstart)
c/7
      call nltest(2,2,1,'rosn','brok',rstart)
c/
      v(29) = amax1(1.0e-7, v(29))
      iv(17) = 30
      iv(18) = 20
c  ***  brown  ***
c/6
c      call nltest(20,4,18,4hbrow,4hn   ,rstart)
c/7
      call nltest(20,4,18,'brow','n   ',rstart)
c/
c
      if (nprob .eq. 0 .or. pu .eq. 0) stop
      call today(datime)
      do 130 k = 1, nprob
         if (mod(k,56) .eq. 1) write(pu, 110) datime, nprob
 110     format(1h1,11x,2a4,2x,2a4,10x,10hsummary of,i4,
     1          22h nl2sol test runs.....,10x,
     2          32h(* = finite-difference jacobian)/
     3          48h0 problem    n   p  niter   nf   ng  iv1  x0scal,5x,
     4          39hfinal f     preldf     nreldf     reldx/)
         j = is(6,k)
         write(pu,120) jtyp(j), name(1,k), name(2,k),
     1                 (is(i,k), i=1,5), irc(k), (rs(i,k), i=1,5)
 120     format(1x,a1,2a4,2i4,i7,2i5,3x,a1,f9.1,e13.3,3e11.3)
 130     continue
c

      stop
c...... last card of nlmain ............................................
      end
