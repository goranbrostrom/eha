      subroutine weibnr(iter, eps, printlevel, nn, ncov, bdim,
     &     time0, time, ind, covar, offset,
     &     beta, loglik, dloglik,
     &     d2loglik, ns, nstra,
     &     conver, fail)

C +++
C     bdim = ncov + 2 * ns
C     beta(bdim + 2 * (j - 1) + 1) = log(lambda(j)), j = 1, ns
C     beta(bdim + 2 * j)           = log(shape(j) ), j = 1, ns
C
C     ncov = No. of covariates (columns in the design matrix = 
C                               rows in 'covar')
C     ns = No. of strata.

      implicit none

      integer iter, printlevel, nn, ncov, bdim
      double precision eps, time0(nn), time(nn), offset(nn)
      double precision covar(ncov, nn), d2loglik(bdim, bdim)
      integer ind(nn), ns, nstra(ns + 1), conver, fail
      double precision beta(bdim)
      double precision loglik, dloglik(bdim)

      integer ione
      parameter (ione = 1)

      double precision one
      parameter (one = 1.d0)

      integer order, i, j, ipfixed

      integer job, itmax, info
      double precision det(2), dnrm2

      double precision db(bdim), L2, pfix

C +++ For dpodi to calculate the inverse only ('det' is a dummy):

      job = ione

C +++ We estimate 'shape', AKA 'p' here
      pfix = 0.d0
      ipfixed = 0
      order = 2

      itmax = iter
      iter = 0

C      call wfunc(order, ipfixed, pfix, bdim, ncov, beta,
C     &     nn, covar, time0, time, ind, offset,
C     &     loglik, dloglik, d2loglik, fail)

      call swfun(order,
     &     bdim, ncov, beta, 
     &     nn, covar, time0, time, ind, 
     &     offset, ns, nstra,
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
C            call intpr("fail in [dpofa]; info:", 22, info, 1) 
            fail = info
            return
         endif
         
C +++ Score test when iter .eq. 0:
C         if (iter .eq. 0) 
C     &        sctest = ddot(ncov, dloglik, ione, db, ione)
            
         L2 = dnrm2(bdim, db, ione)
         if (L2 .lt. eps) conver = 1
         if (printlevel .eq. 1) then
            call intpr("*** Iteration ", 14, iter, 1)
            call dblepr('L2 = ', 5, L2, 1)
            call dblepr("loglik = ", 9, loglik, 1)
         endif

C         if (conver .eq. 0) then
C +++ Update beta:
            call daxpy(bdim, one, db, ione, beta, ione)
            
C            call wfunc(order, ipfixed, pfix, bdim, ncov, beta,
C     &           nn, covar, time0, time, ind, offset,
C     &           loglik, dloglik, d2loglik, fail)

      call swfun(order,
     &     bdim, ncov, beta, 
     &     nn, covar, time0, time, ind, 
     &     offset, ns, nstra,
     &     loglik, dloglik, d2loglik, fail)

      do i = 1, bdim
         dloglik(i) = -dloglik(i)
      enddo
      loglik = -loglik
      
C +++ New iteration
      iter = iter + 1
C         endif

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
C         lambda = exp(beta(bdim))
C         lambda_sd = sqrt(d2loglik(bdim, bdim)) * lambda
         
      else
         fail = info
         return
      endif


      if (printlevel .eq. 1) then
         call intpr('*** Iteration ', 14, iter, 1)
         if (conver .eq. 1) then
            call intpr("Convergence", 11, iter, 0) 
         else
            call intpr("NOTE: No convergence!", 21, iter, 0)
         endif
         call dblepr("loglik = ", 9, loglik, 1)
      endif

      return
      end
