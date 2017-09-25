C Newton - Raphson for 'exponential' regression. Göran Broström (1982-2003).

      subroutine expnr(iter, eps, printlevel, nn, ncov, bdim,
     &     time0, time, ind, covar, offset, pfix,
     &     beta, lambda, lambda_sd, loglik, dloglik, 
     &     d2loglik,
     &     conver, fail)

      implicit none

      integer iter, printlevel, nn, ncov, bdim, conver, fail
      integer ind(nn)
      double precision time0(nn), time(nn), covar(ncov, nn)
      double precision offset(nn), pfix, eps
      double precision beta(bdim), loglik, dloglik(bdim)
      double precision d2loglik(bdim, bdim)
      double precision lambda, lambda_sd

      integer ione
      parameter (ione = 1)

      double precision one
      parameter (one = 1.d0)

      integer order, i, j, ipfixed

      integer job, itmax, info
      double precision det(2), dnrm2

      double precision db(bdim), L2

C +++ For dpodi to calculate the inverse only ('det' is a dummy):

      job = 1

      ipfixed = 1
      order = 2

      itmax = iter
      iter = 0

      call wfunc(order, ipfixed, pfix, bdim, ncov, beta,
     &     nn, covar, time0, time, ind, offset,
     &     loglik, dloglik, d2loglik, fail)
      do i = 1, bdim
         dloglik(i) = -dloglik(i)
      enddo
      loglik = -loglik

      do while ( (iter .lt. itmax) .and. (conver .eq. 0) )

         call dcopy(bdim, dloglik, ione, db, ione) 
         call dpofa(d2loglik, bdim, bdim, info)
         if (info .eq. 0) then
            call dposl(d2loglik, bdim, bdim, db)
         else
C            call intpr('fail in [dpofa]; info = ', -1, info, 1)
            fail = info
            return
         endif
         
C +++ Score test when iter .eq. 0:
C         if (iter .eq. 0) 
C     &        sctest = ddot(ncov, dloglik, ione, db, ione)
            
         L2 = dnrm2(bdim, db, ione)
         if (L2 .lt. eps) conver = 1
         if (printlevel .eq. 1) then
            call intpr(" ", -1, iter, 0)
            call intpr( '*** Iteration ', -1, iter, 1)
            call dblepr( 'L2 = ', -1, L2, 1)
            call dblepr( 'loglik = ', -1, loglik, 1)
         endif

C +++ Update beta:
         call daxpy(bdim, one, db, ione, beta, ione)
            
         call wfunc(order, ipfixed, pfix, bdim, ncov, beta,
     &        nn, covar, time0, time, ind, offset,
     &        loglik, dloglik, d2loglik, fail)
C     wfunc calculates -loglik, -dloglik, -d2loglik
C     -d2loglik is ok (= observed information), but...
         do i = 1, bdim
            dloglik(i) = -dloglik(i)
         enddo
         loglik = -loglik

C +++ New iteration
         iter = iter + 1
      end do

C +++ Done! The afterwork:

C     The inverse of the hessian:
C      call dcopy(ncov, dloglik, ione, db, ione) 
      call dpofa(d2loglik, bdim, bdim, info)
      if (info .eq. 0) then
         call dpodi(d2loglik, bdim, bdim, det, job)
C --- Fill in the 'lower half' of d2loglik:
         do i = 2, bdim
            do j = 1, i - 1
               d2loglik(i, j) = d2loglik(j, i)
            enddo
         enddo
         lambda = exp(beta(bdim))
         lambda_sd = sqrt(d2loglik(bdim, bdim)) * lambda
         
      else
         fail = info
         return
      endif


      if (printlevel .eq. 1) then
         call intpr(" ", 1, iter, 0)
         call intpr( '*** Iteration ', -1, iter, 1)
         if (conver .eq. 1) then
            call intpr( 'Convergence', -1, iter, 0)
         else
            call intpr( 'NOTE: No convergence!', -1, iter, 0)
         endif
         call dblepr('loglik = ', -1, loglik, 1)
      endif

      return
      end
