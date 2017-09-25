C ***
C
C      subroutine gethaz(totrs, ns, 
C     &     antrs, antevents, size,
C     &     totsize, riskset, 
C     &     nn, 
C     &     score, hazard)

      subroutine gethaz(nn, ns, antrs, size, nevents, totsize, 
     &     riskset, score,
     &     totrs, hazard)

      implicit none

      integer nn, ns, totrs, totsize
C      integer antrs(ns), size(ns), nevents(ns)
      integer antrs(ns), size(totrs), nevents(totrs)
      integer riskset(totsize)
      double precision score(nn)
      double precision hazard(totrs) 
C +++ 
C     totevent : Total number of events. [sum(
C     totrs    : Total number of risk sets.
C     ns       : Number of strata.
C
C     antrs     : antrs(i) = No. of risk sets in stratum i, i = 1, ns.
C     antevents : number of events in each riskset.
C     size      : Size of each risk set.
C
C     totsize  : Sum of the risk set sizes.
C     riskset  : pointers to members of risk sets (length totsize).
C
c     nn     : No. of spells.
C +++


C ************************************************************
C     Local:
C
      double precision zero
      parameter (zero = 0.d0) 

      integer rs, j, i, who

      integer indx, rsindx

      double precision sumscore 
C
C *************************************************************

C +++ Calculate score(i), i = 1, nn:

C      call dgemv(trans, nn, antcov, one, covar, nn, beta, ione, zero,  
C     &     score, ione)
C     
C      do i = 1, nn
C         score(i) = exp(score(i))
C      enddo

      rsindx = 0
      indx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
            sumscore = zero
            do i = 1, size(rsindx)
               indx = indx + 1
               who = riskset(indx)
               sumscore = sumscore + score(who)
            enddo
            hazard(rsindx) = nevents(rsindx) / sumscore
         enddo
      enddo
      
      return
      end
