      subroutine parchk(iv, n, nn, p, v)
c
c  ***  check nl2sol (version 2.2) parameters, print changed values  ***
c
      integer iv(1), n, nn, p
      double precision v(1)
c     dimension iv(*), v(*)
c
      external dfault, rmdcon, vcopy
      double precision rmdcon
c dfault -- supplies dfault parameter values.
c rmdcon -- returns machine-dependent constants.
c vcopy  -- copies one vector to another.
c
c  ***  local variables  ***
c
      integer i, iv1, jtolp, k, l, m, nvdflt, pu
c/6
c      real cngd(3), dflt(3), vn(2,27), which(3)
c/7
      character*4 cngd(3), dflt(3), vn(2,27), which(3)
c/
      double precision big, machep, tiny, vk, vm(27), vx(27), zero
c
c  ***  iv and v subscripts  ***
c
      integer dtype, dtype0, d0init, epslon, inits, jtinit, jtol0,
     1        jtol1, oldn, oldnn, oldp, parprt, parsv1, prunit
c
c/6
c      data nvdflt/27/, zero/0.d+0/
c/7
      parameter (nvdflt=27, zero=0.d+0)
c/
c
c/6
c      data dtype/16/, dtype0/29/, d0init/37/, epslon/19/,
c     1     inits/25/, jtinit/39/, jtol0/86/, jtol1/87/,
c     2     oldn/45/, oldnn/46/, oldp/47/, parprt/20/,
c     3     parsv1/51/, prunit/21/
c/7
      parameter (dtype=16, dtype0=29, d0init=37, epslon=19,
     1     inits=25, jtinit=39, jtol0=86, jtol1=87,
     2     oldn=45, oldnn=46, oldp=47, parprt=20,
     3     parsv1=51, prunit=21)
      save big, tiny
c/
c
      data big/0.d+0/, tiny/1.d+0/
c/6
c      data vn(1,1),vn(2,1)/4hepsl,4hon../
c      data vn(1,2),vn(2,2)/4hphmn,4hfc../
c      data vn(1,3),vn(2,3)/4hphmx,4hfc../
c      data vn(1,4),vn(2,4)/4hdecf,4hac../
c      data vn(1,5),vn(2,5)/4hincf,4hac../
c      data vn(1,6),vn(2,6)/4hrdfc,4hmn../
c      data vn(1,7),vn(2,7)/4hrdfc,4hmx../
c      data vn(1,8),vn(2,8)/4htune,4hr1../
c      data vn(1,9),vn(2,9)/4htune,4hr2../
c      data vn(1,10),vn(2,10)/4htune,4hr3../
c      data vn(1,11),vn(2,11)/4htune,4hr4../
c      data vn(1,12),vn(2,12)/4htune,4hr5../
c      data vn(1,13),vn(2,13)/4hafct,4hol../
c      data vn(1,14),vn(2,14)/4hrfct,4hol../
c      data vn(1,15),vn(2,15)/4hxcto,4hl.../
c      data vn(1,16),vn(2,16)/4hxfto,4hl.../
c      data vn(1,17),vn(2,17)/4hlmax,4h0.../
c      data vn(1,18),vn(2,18)/4hdltf,4hdj../
c      data vn(1,19),vn(2,19)/4hd0in,4hit../
c      data vn(1,20),vn(2,20)/4hdini,4ht.../
c      data vn(1,21),vn(2,21)/4hjtin,4hit../
c      data vn(1,22),vn(2,22)/4hdltf,4hdc../
c      data vn(1,23),vn(2,23)/4hdfac,4h..../
c      data vn(1,24),vn(2,24)/4hrlim,4hit../
c      data vn(1,25),vn(2,25)/4hcosm,4hin../
c      data vn(1,26),vn(2,26)/4hdelt,4ha0../
c      data vn(1,27),vn(2,27)/4hfuzz,4h..../
c/7
      data vn(1,1),vn(2,1)/'epsl','on..'/
      data vn(1,2),vn(2,2)/'phmn','fc..'/
      data vn(1,3),vn(2,3)/'phmx','fc..'/
      data vn(1,4),vn(2,4)/'decf','ac..'/
      data vn(1,5),vn(2,5)/'incf','ac..'/
      data vn(1,6),vn(2,6)/'rdfc','mn..'/
      data vn(1,7),vn(2,7)/'rdfc','mx..'/
      data vn(1,8),vn(2,8)/'tune','r1..'/
      data vn(1,9),vn(2,9)/'tune','r2..'/
      data vn(1,10),vn(2,10)/'tune','r3..'/
      data vn(1,11),vn(2,11)/'tune','r4..'/
      data vn(1,12),vn(2,12)/'tune','r5..'/
      data vn(1,13),vn(2,13)/'afct','ol..'/
      data vn(1,14),vn(2,14)/'rfct','ol..'/
      data vn(1,15),vn(2,15)/'xcto','l...'/
      data vn(1,16),vn(2,16)/'xfto','l...'/
      data vn(1,17),vn(2,17)/'lmax','0...'/
      data vn(1,18),vn(2,18)/'dltf','dj..'/
      data vn(1,19),vn(2,19)/'d0in','it..'/
      data vn(1,20),vn(2,20)/'dini','t...'/
      data vn(1,21),vn(2,21)/'jtin','it..'/
      data vn(1,22),vn(2,22)/'dltf','dc..'/
      data vn(1,23),vn(2,23)/'dfac','....'/
      data vn(1,24),vn(2,24)/'rlim','it..'/
      data vn(1,25),vn(2,25)/'cosm','in..'/
      data vn(1,26),vn(2,26)/'delt','a0..'/
      data vn(1,27),vn(2,27)/'fuzz','....'/
c/
c
      data vm(1)/1.0d-3/, vm(2)/-0.99d+0/, vm(3)/1.0d-3/, vm(4)/1.0d-2/,
     1     vm(5)/1.2d+0/, vm(6)/1.d-2/, vm(7)/1.2d+0/, vm(8)/0.d+0/,
     2     vm(9)/0.d+0/, vm(10)/1.d-3/, vm(11)/-1.d+0/, vm(15)/0.d+0/,
     3     vm(16)/0.d+0/, vm(19)/0.d+0/, vm(20)/-10.d+0/, vm(21)/0.d+0/,
     4     vm(23)/0.d+0/, vm(24)/1.d+10/, vm(27)/1.01d+0/
      data vx(1)/0.9d+0/, vx(2)/-1.d-3/, vx(3)/1.d+1/, vx(4)/0.8d+0/,
     1     vx(5)/1.d+2/, vx(6)/0.8d+0/, vx(7)/1.d+2/, vx(8)/0.5d+0/,
     2     vx(9)/0.5d+0/, vx(10)/1.d+0/, vx(11)/1.d+0/, vx(14)/0.1d+0/,
     3     vx(15)/1.d+0/, vx(16)/1.d+0/, vx(18)/1.d+0/, vx(22)/1.d+0/,
     4     vx(23)/1.d+0/, vx(25)/1.d+0/, vx(26)/1.d+0/, vx(27)/1.d+2/
c
c/6
c      data cngd(1),cngd(2),cngd(3)/4h---c,4hhang,4hed v/,
c     1     dflt(1),dflt(2),dflt(3)/4hnond,4hefau,4hlt v/
c/7
      data cngd(1),cngd(2),cngd(3)/'---c','hang','ed v'/,
     1     dflt(1),dflt(2),dflt(3)/'nond','efau','lt v'/
c/
c
c.......................................................................
c
      if (iv(1) .eq. 0) call dfault(iv, v)
      pu = iv(prunit)
      iv1 = iv(1)
      if (iv1 .ne. 12) go to 30
         if (nn .ge. n .and. n .ge. p .and. p .ge. 1) go to 20
              iv(1) = 16
              if (pu .ne. 0) write(pu,10) nn, n, p
 10           format(30h0///// bad nn, n, or p... nn =,i5,5h, n =,i5,
     1               5h, p =,i5)
              go to 999
 20      k = iv(21)
         call dfault(iv(21), v(33))
         iv(21) = k
         iv(dtype0) = iv(dtype+20)
         iv(oldn) = n
         iv(oldnn) = nn
         iv(oldp) = p
         which(1) = dflt(1)
         which(2) = dflt(2)
         which(3) = dflt(3)
         go to 80
 30   if (n .eq. iv(oldn) .and. nn .eq. iv(oldnn) .and. p .eq. iv(oldp))
     1                       go to 50
         iv(1) = 17
         if (pu .ne. 0) write(pu,40) iv(oldnn), iv(oldn), iv(oldp), nn,
     1                               n, p
 40      format(30h0///// (nn,n,p) changed from (,i5,1h,,i5,1h,,i3,
     1          6h) to (,i5,1h,,i5,1h,,i3,2h).)
         go to 999
c
 50   if (iv1 .le. 11 .and. iv1 .ge. 1) go to 70
         iv(1) = 50
         if (pu .ne. 0) write(pu,60) iv1
 60      format(15h0/////  iv(1) =,i5,28h should be between 0 and 12.)
         go to 999
c
 70   which(1) = cngd(1)
      which(2) = cngd(2)
      which(3) = cngd(3)
c
 80   if (big .gt. tiny) go to 90
         tiny = rmdcon(1)
         machep = rmdcon(3)
         big = rmdcon(6)
         vm(12) = machep
         vx(12) = big
         vm(13) = tiny
         vx(13) = big
         vm(14) = machep
         vm(17) = tiny
         vx(17) = big
         vm(18) = machep
         vx(19) = big
         vx(20) = big
         vx(21) = big
         vm(22) = machep
         vx(24) = rmdcon(5)
         vm(25) = machep
         vm(26) = machep
 90   m = 0
      if (iv(inits) .ge. 0 .and. iv(inits) .le. 2) go to 110
         m = 18
         if (pu .ne. 0) write(pu,100) iv(inits)
 100     format(25h0/////  inits... iv(25) =,i4,20h should be between 0,
     1          7h and 2.)
 110  k = epslon
      do 140 i = 1, nvdflt
         vk = v(k)
         if (vk .ge. vm(i) .and. vk .le. vx(i)) go to 130
              m = k
              if (pu .ne. 0) write(pu,120) vn(1,i), vn(2,i), k, vk,
     1                                    vm(i), vx(i)
 120          format(8h0/////  ,2a4,5h.. v(,i2,3h) =,d11.3,7h should,
     1               11h be between,d11.3,4h and,d11.3)
 130     k = k + 1
 140     continue
c
      if (iv1 .eq. 12 .and. v(jtinit) .gt. zero) go to 170
c
c  ***  check jtol values  ***
c
      jtolp = jtol0 + p
      do 160 i = jtol1, jtolp
         if (v(i) .gt. zero) go to 160
         k = i - jtol0
         if (pu .ne. 0) write(pu,150) k, i, v(i)
 150     format(12h0///// jtol(,i3,6h) = v(,i3,3h) =,d11.3,
     1          20h should be positive.)
         m = i
 160     continue
c
 170  if (m .eq. 0) go to 180
         iv(1) = m
         go to 999
c
 180  if (pu .eq. 0 .or. iv(parprt) .eq. 0) go to 999
      if (iv1 .ne. 12 .or. iv(inits) .eq. 0) go to 200
         m = 1
         write(pu,190) iv(inits)
 190     format(22h0nondefault values..../20h inits..... iv(25) =,i3)
 200  if (iv(dtype) .eq. iv(dtype0)) go to  210
         if (m .eq. 0) write(pu,215) which
         m = 1
         write(pu,205) iv(dtype)
 205     format(20h dtype..... iv(16) =,i3)
 210  k = epslon
      l = parsv1
      do 240 i = 1, nvdflt
         if (v(k) .eq. v(l)) go to 230
              if (m .eq. 0) write(pu,215) which
 215          format(1h0,3a4,9halues..../)
              m = 1
              write(pu,220) vn(1,i), vn(2,i), k, v(k)
 220          format(1x,2a4,5h.. v(,i2,3h) =,d15.7)
 230     k = k + 1
         l = l + 1
 240     continue
      iv(dtype0) = iv(dtype)
      call vcopy(nvdflt, v(parsv1), v(epslon))
      if (iv1 .ne. 12) go to 999
         if (v(jtinit) .gt. zero) go to 260
              jtolp = jtol0 + p
              write(pu,250) (v(i), i = jtol1, jtolp)
 250          format(24h0(initial) jtol array.../(1x,6d12.3))
 260     if (v(d0init) .gt. zero) go to 999
              k = jtol1 + p
              l = k + p - 1
              write(pu,270) (v(i), i = k, l)
 270          format(22h0(initial) d0 array.../1x,6d12.3)
c
 999  return
c  ***  last card of parchk follows  ***
      end
