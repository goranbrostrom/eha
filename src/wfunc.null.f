C ***
C
      subroutine wfuncnull(order, ipfixed, pfix,
     &     bdim, b, nn, time0, time, ind,
     &     f, fp, fpp, iok)

C ******************************************************************** C
C
C     Calculates f, the negative of a Weibull loglihood function without
C     covariates. Eventually calculates the first (fp) and second (fpp)
C     order derivatives as well; see the input variable 'order' below.
C
C     Note: No stratification; or only calculations in one stratum.
C ******************************************************************** C

      implicit none

C      include 'param.inc'
C      include 'optcom.inc'

      integer order, ipfixed, bdim, nn, iok
   
      double precision pfix
      double precision time0(nn), time(nn)
      integer ind(nn)

      double precision b(bdim), f, fp(bdim)
      double precision fpp(bdim, bdim)

C ***
C     INPUT:
C
C     order:  negative: return
C             0       : calculate only f.
C             1       : calculate  f and fp.
C             >= 2    : calculate f, fp and fpp.
C     b:    b(1) = log(lambda); b(2) = log(p)
C
C     Output:
C
C     f     : negative of log likelihood function,
C     fp    : negative of score vector.
C     fpp   : negative of hessian.
C
C     f, fp and fpp are calculated at lambda (and p).
C ***
C
C ***
C     Baseline survival function = exp(-(lambda*time)**p).
C     lambda = b(1) = exp(alfa), p = b(2) = exp(gamma).
C ***
      double precision d, dy

      double precision s, sy, syy

      double precision alfa, gamma, p, ap, pps, ppsy
      integer j1, id, index

      logical pfixed, ord1, ord2, ok

      double precision zero, one, oflimit
      parameter (zero = 0.d0, one = 1.d0, oflimit = 500.d0)

C ***
C     If order is negative, do nothing!
C ***

C     To keep the compiler happy (iok is not used here):
      iok = 0
C
      if (order .lt. 0) return

      pfixed = (ipfixed .ne. 0)

      ord1 = (order .ge. 1)
      ord2 = (order .ge. 2)
C +++
C     Calculate the "help sums", i.e. s, sy, etc.
C +++
      if (pfixed) then
         ok = .true.
         p = pfix
         gamma = log(p)
      else
         gamma = b(2)

C +++
C     Protect against overflow and underflow in the EXP function.
C     Skip that (030525).
C +++
C         ok = (abs(gamma) .le. oflimit)
C         if (.not. ok) return
C         if (.not. ok) then
C            call intpr('order = ', 8, order, 1)
C            call dblepr('gamma = ', 8, gamma, 1)
C            call dblepr('p = ', 4, p, 1)
C            call rexit(" abs(gamma) > oflimit!!")
C         endif
C
         p = exp(gamma)
      endif

      alfa = b(1)
      ap = alfa * p
C +++
C     Calculate d, dy, dz if first time with this data set.
C
C     d:      total number of observed deaths.
C     dy:     sum of d(i) * log[time(i)] over all individuals.
C +++
C      if (iter .le. 1) then
      id = 0
      d = zero
      dy = zero

      do 7 j1 = 1, nn
         if ( (ind(j1) .eq. 1) .and. (time(j1) .gt. zero) )then
            id = id + 1
            dy = dy + log(time(j1))
         endif
 7    continue
      d = dble(id)

      call getsums_null(ord1, ord2, alfa, p, pfixed, 
     +        nn, time, time0, ind,
     +        s, sy, syy)

C ***
C     Calculate f.
C ***
      f = dy * (one - p) + s - d * (gamma + ap)

C ***
C     If order > 0, calculate first derivatives.
C ***
      if (ord1)  then
       fp(1) = p * (s - d)
       if (.not. pfixed) fp(2) = p * (alfa * s + sy) - d * (one +
     +            ap) - dy * p


C ***
C     If order > 1, calculate second derivatives.
C ***
       if (ord2)  then
          index = 0

          index = index + 1
          pps = p * s * p
          fpp(1, 1) = pps

            if (.not. pfixed) then
               ppsy = p * sy * p
               index = index + 1
               fpp(1, 2) = alfa * pps + ppsy + fp(1)
               index = index + 1
               fpp(2, 2) = alfa * (alfa * pps + 2.d0 * ppsy) + 
     $                      p * syy * p + d + fp(2)
C +++ Fill in missing triangle:
               fpp(2, 1) = fpp(1, 2)
            endif
         endif
C +++    order .ge. 2
      endif
C +++ order .ge. 1

      return
      end
