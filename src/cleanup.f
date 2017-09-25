      subroutine cleanup(ocovar, oenter, oexit, oevent, oid,
     1     ncov, onrec, onn, eps,
     2     rec, covar, enter, exi, event, id)

C     ocovar = old covariate matrix
C     oenter = old enter
C     oexit = old exit
C     oevent = old event
C     oid = old id
C     ncov = number of covariates (old and new)
c     onrec = # of old records
C     onn = # of old individuals ('length(unique(id))')
C     rec = # of new records
C     covar = new covariate matrix
C     enter = new enter
C     exi = new exit
C     event = new event
C     id = new id

      implicit none

      integer ncov
      integer onrec, onn
      integer rec

      double precision ocovar(ncov, onrec), covar(ncov, onrec)
      double precision oenter(onrec), oexit(onrec)
      double precision enter(onrec), exi(onrec), eps
      integer oevent(onrec), event(onrec), oid(onrec), id(onrec)

      integer i, pers, oldpers, startrec, start

C      logical equal, cequal

      integer antrec(onn)

      rec = 0
      startrec = 0

      oldpers = oid(1)

      do i = 1, onn
         antrec(i) = 0
      enddo

      do i = 1, onrec
         pers = oid(i)
         antrec(pers) = antrec(pers) + 1
      enddo

      start = 1
      call persout(oid(start), oenter, oexit, oevent, ncov, antrec(1), 
     &     ocovar(1, start), onrec, id, enter, exi, event, covar, rec, 
     &     eps)
 
      
      do i = 2, onn
         start = start + antrec(i - 1)
         call persout(oid(start), oenter(start), oexit(start), 
     &        oevent(start), 
     &        ncov, antrec(i), ocovar(1, start), onrec, id, 
     &        enter, exi, event, covar, rec, eps)
 
      enddo

      return

      end

C ***
C
      logical function equal(x, y, eps)

      implicit none

      double precision x, y
      double precision  eps

      equal = (abs(x - y) .lt. eps)
      
      return

      end

C ***
C
      logical function cequal(n, x, y, eps)

      implicit none

      integer n

      double precision x(n), y(n), eps
      integer i
      logical res

      res = .TRUE.
      i = 0
      do while (res .and. (i .lt. n))
         i = i + 1
         res = (abs(x(i) - y(i)) .lt. eps)
      enddo
      
      cequal = res

      return

      end
C ***
C
      subroutine putrec(rec, pers, id, oenter, enter, oexit, exi,
     1                    oevent, event, ocovar, covar, ncov, onrec)

      implicit none

      integer rec, pers, ncov, onrec
      double precision oenter, oexit, ocovar(ncov)
      integer oevent
      double precision enter(onrec), exi(onrec), covar(ncov, onrec)
      integer id(onrec), event(onrec)

      integer i

      enter(rec) = oenter
      exi(rec) = oexit
      event(rec) = oevent
      id(rec) = pers

      do i = 1, ncov
         covar(i, rec) = ocovar(i)
      enddo

      return

      end

C ***
C
      subroutine persout(oid, oenter, oexit, oevent, ncov, dim, ocovar, 
     1     onrec, id, enter, exi, event, covar, rec, eps)

      implicit none

      integer oid, ncov, dim, onrec, rec
      integer oevent(dim), event(onrec), id(onrec) 
      double precision oenter(dim), oexit(dim), enter(onrec), exi(onrec)
      double precision ocovar(ncov, dim), covar(ncov, onrec), eps

      integer i, antdead

      logical equal, cequal

      rec = rec + 1

      if (oevent(1) .gt. 0) then
         antdead = 1
      else
         antdead = 0
      endif

      call putrec(rec, oid, id, oenter(1), enter, 
     &     oexit(1), exi, oevent(1), event, 
     &     ocovar(1, 1), covar, ncov, onrec)

      if (antdead .eq. 1) return

      do i = 2, dim

         if (oevent(i) .gt. 0) then
            antdead = 1
         else
            antdead = 0
         endif
         
         if (equal( exi(rec), oenter(i), eps )) then
            if (cequal(ncov, covar(1, rec), ocovar(1, i), eps )) then
               exi(rec) = oexit(i)
               event(rec) = oevent(i)
            else
               rec = rec + 1
               call putrec(rec, oid, id, oenter(i), enter, 
     &              oexit(i), exi, oevent(i), event, 
     &              ocovar(1, i), covar, ncov, onrec)
            endif
         elseif (exi(rec) .lt. enter(i)) then
            rec = rec + 1
            call putrec(rec, oid, id, oenter(i), enter, 
     &           oexit(i), exi, oevent(i), event, 
     &           ocovar(1, i), covar, ncov, onrec)
         else
C     *** This is the critical case; overlapping spells?
            if ( (oexit(i) .ge. exi(rec)) .or. (antdead .eq. 1) ) then
               exi(rec) = oenter(i)
               if (cequal(ncov, covar(1, rec), ocovar(1, i), eps )) then
                  exi(rec) = oexit(i)
                  event(rec) = oevent(i)
               else
                  rec = rec + 1
                  call putrec(rec, oid, id, oenter(i), enter, 
     &                 oexit(i), exi, oevent(i), event, 
     &                 ocovar(1, i), covar, ncov, onrec)
               endif
               
            endif
         endif
         if (antdead .eq. 1) return
      enddo

      return

      end
