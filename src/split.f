      subroutine split(xx, nn, colx, yy, nny, coly, 
     &     nn_out, ind_iv, cuts, n_years)

      implicit none

      integer nn, colx, nny, coly, n_years
      integer nn_out(nn), ind_iv(nn, 2)
      double precision xx(nn, colx), yy(nny, coly)
      double precision cuts(n_years + 1)

      integer i, j, begrow, endrow, row
      double precision iv_length

      iv_length = cuts(2) - cuts(1)

      begrow = 0
      endrow = 0
      do j = 1, nn
         begrow = endrow + 1
         endrow = endrow + nn_out(j)
         if (begrow .eq. endrow) then
            do i = 1, colx
               yy(begrow, i) = xx(j, i)
            enddo
            yy(begrow, colx + 1) = ind_iv(j, 1)
         else
            do row = begrow, endrow
               yy(row, 4) = xx(j, 4)
               yy(row, 5) = xx(j, 5)
            enddo

            row = begrow
            yy(row, 1) = xx(j, 1)
            yy(row, 2) = cuts(ind_iv(j, 1) + 1) - xx(j, 4)
            yy(row, 3) = 0.0
            yy(row, 6) = ind_iv(j, 1)

            if ( (endrow - begrow) .ge. 2) then
               do row = (begrow + 1), (endrow - 1)
                  yy(row, 1) = yy(row - 1, 2)
                  yy(row, 2) = yy(row - 1, 2) + iv_length
                  yy(row, 6) = ind_iv(j, 1) + row - begrow
               enddo
            endif

            row = endrow
            yy(row, 1) = cuts(ind_iv(j, 2)) - xx(j, 4)
            yy(row, 2) = xx(j, 2)
            yy(row, 3) = xx(j, 3)
            yy(row, 6) = ind_iv(j, 2)

         endif
      enddo

      return

      end
         
