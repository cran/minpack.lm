      double precision function ss(n,x)
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
      data zero /0.0d0/
      ss = zero
      do 10 i = 1, n
         ss = ss + (x(i))**2
   10   continue
      return
c     
c     last card of function ss.
c
      end
