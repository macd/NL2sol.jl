      double precision function rmdcon(k)
c
c  ***  return machine dependent constants used by nl2sol  ***
c
c +++  comments below contain data statements for various machines.  +++
c +++  to convert to another machine, place a c in column 1 of the   +++
c +++  data statement line(s) that correspond to the current machine +++
c +++  and remove the c from column 1 of the data statement line(s)  +++
c +++  that correspond to the new machine.                           +++
c
      integer k
c
c  ***  the constant returned depends on k...
c
c  ***        k = 1... smallest pos. eta such that -eta exists.
c  ***        k = 2... square root of 1.001*eta.
c  ***        k = 3... unit roundoff = smallest pos. no. machep such
c  ***                 that 1 + machep .gt. 1 .and. 1 - machep .lt. 1.
c  ***        k = 4... square root of 0.999*machep.
c  ***        k = 5... square root of 0.999*big (see k = 6).
c  ***        k = 6... largest machine no. big such that -big exists.
c
      double precision big, eta, machep
c/+
      double precision dsqrt
c/
      double precision one001, pt999
c
      data one001/1.001d0/, pt999/0.999d0/
c
c  +++  ibm 360, ibm 370, or xerox  +++
c
c      data big/z7fffffffffffffff/, eta/z0010000000000000/,
c     1     machep/z3410000000000000/
c
c  +++  data general  +++
c
c     data big/0.7237005577d+76/, eta/0.5397605347d-78/,
c    1     machep/2.22044605d-16/
c
c  +++  dec 11  +++
c
c     data big/1.7d+38/, eta/2.938735878d-39/, machep/2.775557562d-17/
c
c  +++  hp3000  +++
c
c     data big/1.157920892d+77/, eta/8.636168556d-78/,
c    1     machep/5.551115124d-17/
c
c  +++  honeywell  +++
c
c     data big/1.69d+38/, eta/5.9d-39/, machep/2.1680435d-19/
c
c  +++  dec10  +++
c
c     data big/"377777100000000000000000/,
c    1     eta/"002400400000000000000000/,
c    2     machep/"104400000000000000000000/
c
c  +++  burroughs  +++
c
c     data big/o0777777777777777,o7777777777777777/,
c    1     eta/o1771000000000000,o7770000000000000/,
c    2     machep/o1451000000000000,o0000000000000000/
c
c  +++  control data  +++
c
c
c     data big/37767777777777777777b,37167777777777777777b/,
c    1     eta/00014000000000000000b,00000000000000000000b/,
c    2     machep/15614000000000000000b,15010000000000000000b/
c
c  +++  prime  +++
c
c     data big/1.0d+9786/, eta/1.0d-9860/, machep/1.4210855d-14/
c
c  +++  univac  +++
c
c     data big/8.988d+307/, eta/1.2d-308/, machep/1.734723476d-18/
c
c  +++  vax  +++
c
c     data big/1.7d+38/, eta/2.939d-39/, machep/1.3877788d-17/
C
c  +++  Modern day x64 Linux  +++
c
      parameter (machep=2.220446049250313D-16)
      parameter (eta=2.2250738585072014D-308)
      parameter (big=1.7976931348623157D+308) 
c
c-------------------------------  body  --------------------------------
c
      go to (10, 20, 30, 40, 50, 60), k
c
 10   rmdcon = eta
      go to 999
c
 20   rmdcon = dsqrt(one001*eta)
      go to 999
c
 30   rmdcon = machep
      go to 999
c
 40   rmdcon = dsqrt(pt999*machep)
      go to 999
c
 50   rmdcon = dsqrt(pt999*big)
      go to 999
c
 60   rmdcon = big
c
 999  return
c  ***  last card of rmdcon follows  ***
      end
