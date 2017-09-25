C ***
C
      subroutine hazards(totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, riskset, 
     &     nn, 
C        covar, strata,
     &     score, hazard)

C +++ 
C     totevent : Total number of events.
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
C     antcov : No. of covariates.
C     covar  : matrix of covariates (nn x antcov).
C     strata : Vector of stratum indicators.
C     beta   : Vector of coefficients.
C
C     score, sumdscore, sumd2score: 'Work areas', avoiding local
C                                   dynamic memory allocation.
C +++

      implicit none

      integer totrs, ns, totsize, nn
      integer antrs(ns), antevents(totrs), size(totrs)
      integer riskset(totsize)
C      double precision covar(nn, antcov)
C      integer strata(nn)

      double precision hazard(totrs)

C +++ Work areas:
      double precision score(nn)
C ************************************************************
C     Local:
C
      double precision zero, one
      parameter (zero = 0.d0, one = 1.d0) 
      integer ione
      parameter (ione = 1)
      character*1 trans
      parameter (trans = 'N')

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
            sumscore = 0.d0
            do i = 1, size(rsindx)
               indx = indx + 1
               who = riskset(indx)
               sumscore = sumscore + score(who)
            enddo
            hazard(rsindx) = antevents(rsindx) / sumscore
         enddo
      enddo
      
      return
      end
