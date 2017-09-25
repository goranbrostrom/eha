C ***
C
      subroutine martres(totrs, ns, 
     &     antrs, antevents, size,
     &     totsize, riskset, 
     &     nn,
     &     score, hazard, resid)

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
C     hazard : Nelson-Aalen 'atoms'.
C     beta   : Output: martingale residuals.
C +++

      implicit none

      integer totrs, ns, totsize, nn
      integer antrs(ns), antevents(totrs), size(totrs)
      integer riskset(totsize)

      double precision score(nn), hazard(totrs), resid(nn)

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

      double precision haz
C
C *************************************************************

      do i = 1, nn
         resid(i) = 0.d0
      enddo

      rsindx = 0
      indx = 0
      do rs = 1, ns
         do j = 1, antrs(rs)
            rsindx = rsindx + 1
            haz = hazard(rsindx)
            do i = 1, antevents(rsindx)
               indx = indx + 1
               who = riskset(indx)
               resid(who) = resid(who) + 1.d0 - haz * score(who)
            enddo
            do i = antevents(rsindx) + 1, size(rsindx)
               indx = indx + 1
               who = riskset(indx)
               resid(who) = resid(who) - haz * score(who)
            enddo
         enddo
      enddo

      return
      end
