C ***
C
      subroutine wfunc(order, ipfixed, pfix, bdim, mb, b, 
     &     nn, z, time0, time, ind, offset,
     &     f, fp, fpp, iok)

C ******************************************************************** C
C
C     Calculates f, the negative of a Weibull loglihood function with
C     covariates. Eventually calculates the first (fp) and second (fpp)
C     order derivatives as well; see the input variable 'order' below.
C
C     Note: No stratification; or only calculations in one stratum.
C ******************************************************************** C

      implicit none

C      include 'param.inc'
C      include 'optcom.inc'

      integer order, ipfixed,bdim, mb, nn, iok
   
      double precision pfix
      double precision z(mb, nn), time0(nn), time(nn), offset(nn)
      integer ind(nn)
      double precision b(bdim)
      double precision f, fp(bdim)
      double precision fpp(bdim, bdim)

C ***
C     INPUT:
C
C     order:  negative: return
C             0       : calculate only f.
C             1       : calculate  f and fp.
C             >= 2    : calculate f, fp and fpp.
C     bdim  : dimension of b, including ln(lambda),
c             and ln(p), if .not.pfixed.
C     b     : vector of beta-estimates,
C             b(k+1) = alfa = ln(lambda), b(k+2) = gamma = ln(p).
C
C     In COMMON BLOCK funccom.inc:
C
C     nn:     No. of records, ind = 0, 1, 2.
C
C     pfixed: .true., if p is regarded as a fixed constant.
C             .false., if we maximize w.r.t. p also.
C
C
C
C     Output:
C
C     f     : negative of log likelihood function,
C     fp    : negative of score vector.
C     fpp   : negative of hessian.
C
C     f, fp and fpp are calculated at b (and p).
C ***
C
C ***
C     Baseline survival function = exp(-(lambda*time)**p).
C     lambda = b(k+1) = exp(alfa), p = b(k+2) = exp(gamma).
C ***
      double precision d, dy, dz(mb)
C      save d, dy, dz

      double precision s, sy, syy, sz(mb), syz(mb)
      double precision szz(mb * (mb + 1) / 2)
      double precision alfa, gamma, bdz, p, ap, pps, ppsy
      integer k, kp1, kp2, j1, j2, id, index

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
         k = bdim - 1
         p = pfix
         gamma = log(p)
      else
         k = bdim - 2
         gamma = b(bdim)

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

      kp1 = k + 1
      kp2 = kp1 + 1
      alfa = b(kp1)
      ap = alfa * p
C +++
C     Calculate d, dy, dz if first time with this data set.
C
C     d:      total number of observed deaths.
C     dy:     sum of d(i) * log[time(i)] over all individuals.
C     dz:     dz(j)  is sum of z(j, i) * d(i) over
C            i (individuals), j = 1,...k.
C +++
C      if (iter .le. 1) then
      id = 0
      d = zero
      dy = zero
      call dcopy(k, zero, 0, dz, 1)
C     do 3 j1 = 1, k
C     dz(j1) = zero
C     3    continue
      do 7 j1 = 1, nn
C     call GetRec(j1, wtime, th0,
C     $              wz, woffset, wcommun, wind, wstratum, wrank, ok)
C            if (.not. ok) stop ' Error in getting record!!'
         if ( (ind(j1) .eq. 1) .and. (time(j1) .gt. zero) )then
C          if ((wind .eq. 1) .and. (wtime .gt. zero)) then
            id = id + 1
            dy = dy + log(time(j1))
            
            do 6 j2 = 1, k
               dz(j2) = dz(j2) + z(j2, j1)
 6          continue
         endif
 7    continue
      d = dble(id)
C      endif

      call getsums(ord1, ord2, mb, k, b, alfa, p, pfixed, 
     +        nn, time, time0, ind, z, offset,
     +        s, sy, syy, sz, syz, szz)
      bdz = zero

      do 10 j1 = 1, k
         bdz = bdz + b(j1) * dz(j1)
   10 continue
C ***
C     Calculate f.
C ***
      f = dy * (one - p) + s - d * (gamma + ap) - bdz

C ***
C     If order > 0, calculate first derivatives.
C ***
      if (ord1)  then
       do 20 j1 = 1, k
          fp(j1) = sz(j1) - dz(j1)
   20       continue
       fp(kp1) = p * (s - d)
       if (.not. pfixed) fp(kp2) = p * (alfa * s + sy) - d * (one +
     +            ap) - dy * p


C ***
C     If order > 1, calculate second derivatives.
C ***
       if (ord2)  then
          index = 0

          do 40 j2 = 1, k
             do 30 j1 = 1, j2
              index = index + 1
C              fpp(index) = szz(index)
              fpp(j1, j2) = szz(index)
   30             continue
   40          continue

            do 50 j1 = 1, k
             index = index + 1
C             fpp(index) = p * sz(j1)
             fpp(j1, k+1) = p * sz(j1)
   50       continue

          index = index + 1
          pps = p * s * p
          fpp(k + 1, k + 1) = pps

            if (.not. pfixed) then
               do 60 j1 = 1, k
                  index = index + 1
                  fpp(j1, k + 2) = p * (alfa * sz(j1) + syz(j1))
   60          continue

               ppsy = p * sy * p
               index = index + 1
               fpp(k + 1, k + 2) = alfa * pps + ppsy + fp(kp1)
               index = index + 1
               fpp(k + 2, k + 2) = alfa * (alfa * pps + 2.d0 * ppsy) + 
     $                      p * syy * p + d + fp(kp2)
            endif
C +++ Fill in missing triangle:
            do j1 = 2, bdim
               do j2 = 1, j1 - 1
                  fpp(j1, j2) = fpp(j2, j1)
               enddo
            enddo
         endif
C +++    order .ge. 2
      endif
C +++ order .ge. 1

      return
      end
