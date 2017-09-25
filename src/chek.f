      subroutine chek(n, n_ind, id_size,
     &     enter, exi, event, eps, sane)

      implicit none

      integer n, n_ind, id_size(n_ind), event(n)
      double precision enter(n), exi(n), eps

      integer sane(n_ind)

      logical ok
      integer i, start

      start = 1
      do i = 1, n_ind
         call check_id(id_size(i), 
     &        enter(start), exi(start), event(start), eps, ok)
         if (ok) then
            sane(i) = 1
         else
            sane(i) = 0
         endif
         start = start + id_size(i)
      enddo

      return
      end

      subroutine check_id(n, enter, exi, event, eps, ok)

C +++ Here, n = No. of records for _this_ individual.

      implicit none

      integer n, event(n)
      double precision enter(n), exi(n), eps
      logical ok

      integer i, tot_ev

      ok = ((exi(1) - enter(1)) .ge. eps)
      if (n .eq. 1) then
         ok = ok .and. ((event(1) .eq. 0) .or. (event(1) .eq. 1))
      else
         i = 2
         tot_ev = event(1)
         do while (ok .and. (i .lt. n))
            ok = ((exi(i) - enter(i)) .ge. eps)
            ok = ok .and. (event(i) .eq. 0)
            ok = ok .and. (exi(i-1) .le. enter(i)) 
            i = i + 1
         enddo
         if (ok) then
            ok = ((exi(n) - enter(n)) .ge. eps)
            ok = ok .and. ((event(n) .eq. 0) .or. (event(n) .eq. 1))
            ok = ok .and. (exi(n-1) .le. enter(n)) 
         endif
      endif
            
      return
      end
