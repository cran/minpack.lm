      double precision function getsse(n,x)
      integer n
      double precision x(n)
c     **********
c
c     function ss
c
c     given an n-vector x, this function calculates the
c     sum of squares of x.
c
c     the function statement is
c
c       double precision function ss(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c     k. mullen, march 2008
c
c     **********
      integer i
      double precision s1
      data zero /0.0d0/
      s1 = zero
      do 90 i = 1, n
         s1 = s1 + x(i)**2
   90    continue
      ss = s1 
  130 continue
      return
c
c     last card of function ss.
c
      end
