C ***
C
      subroutine update_sums(ord1, ord2, k, wind, wtime, wz, 
     &     woffset, pfixed, p, alfa, b, 
     &     s, sy, syy, sz, syz, szz)


      implicit none

      logical ord1, ord2
      integer k, wind
      double precision b(k), woffset
      logical pfixed
      double precision p, alfa
      double precision wtime, wz(k), s, sy, syy, sz(k), syz(k)
      double precision szz(k * (k + 1) / 2)

C +++ Local:
      double precision zero
      parameter (zero = 0.d0)
      double precision y, zma(k), zb, ddot, epy
      double precision tmp1, tmp2
      integer index, j2

      call dcopy(k, wz, 1, zma, 1)
      zb = ddot(k, b, 1, zma, 1)
      zb = zb + woffset

C     do 30 j1 = 1, k
C     zma(j1) = z(j1, i)
C     zb = zb + b(j1) * zma(j1)
C     30       continue
      
      y = log(wtime)
      
      epy = exp(p * (y + alfa) + zb)
C      if (epy .lt. -250.d0) then
C         epy = zero
C      elseif (epy .le. 250.d0) then
C         epy = exp(epy)
C      else
C         epy = exp(250.d0)
C      endif
      
C     FY:          if (ind(i) .eq. 2) epy = -epy
      if (wind .eq. 2) epy = -epy
              
      s = s + epy
              
      if (ord1) then
         call daxpy(k, epy, zma, 1, sz, 1)
C     do 40 j1 = 1, k
C     sz(j1) = sz(j1) + zma(j1) * epy
C     40          continue
                 
         if (.not. pfixed) then
            tmp1 = y * epy
            sy = sy + tmp1
         endif
         
         if (ord2) then
            index = 1
            do 60 j2 = 1, k
               tmp2 = epy * zma(j2)
               call daxpy(j2, tmp2, zma, 1, szz(index), 1)
               index = index + j2
C     do 50 j1 = 1, j2
C     index = index + 1
C     szz(index) = szz(index) + zma(j1)*zma(j2)*epy
C     50                continue
 60         continue
                    
            if (.not. pfixed) then
               call daxpy(k, tmp1, zma, 1, syz, 1)
C     do 70 j1 = 1, k
C     syz(j1) = syz(j1) + zma(j1) * tmp1
C     70                continue
               syy = syy + tmp1 * y
            endif
         endif
c     +++          order .ge. 2
      endif
c     +++       order .ge. 1
              
      return
      end

C ***
C
      subroutine getsums(ord1, ord2, mb, k, b, alfa, p, pfixed,
     +     nn, time, time0, ind, z, offset,
     +     s, sy, syy, sz, syz, szz)

      implicit none

C      include 'param.inc'
C      include 'optcom.inc'

      logical ord1, ord2, pfixed

      integer mb, k, nn, ind(nn)
      double precision time(nn), time0(nn), z(mb, nn), offset(nn)
      double precision b(k), alfa, p, s, sy, syy, sz(k), syz(k)
      double precision szz(k * (k + 1) / 2)

C +++ Local:
      integer i, index

      double precision wtime
      integer wind

      double precision zero
      parameter (zero = 0.d0)

      s = zero
      sy = zero
      syy = zero

C      do 10 j1 = 1, k
C       sz(j1) = zero
C       syz(j1) = zero
C   10 continue
      call dcopy(k, zero, 0, sz, 1)
      call dcopy(k, zero, 0, syz, 1)

      index = k * (k + 1) / 2 
C      do 20 j1 = 1, index
C       szz(j1) = zero
C   20 continue
      call dcopy(index, zero, 0, szz, 1)

      do 100 i = 1, nn
C         call GetRec(i, wtime, th0,
C     $              wz, woffset, wcommun, wind, wstratum, wrank, ok)
         wtime = time0(i)
         if  (wtime .gt. zero) then
            wind = 2
            call update_sums(ord1, ord2, k, wind, wtime, z(1, i),
     &           offset(i), pfixed, p, alfa, b,
     &           s, sy, syy, sz, syz, szz)
         endif
         wtime = time(i)
         wind = ind(i)
         call update_sums(ord1, ord2, k, wind, wtime, z(1, i),
     &        offset(i), pfixed, p, alfa, b, 
     &        s, sy, syy, sz, syz, szz)
            
  100 continue

      return
      end
