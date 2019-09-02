        SUBROUTINE ghq(n, x, w, modified)
C
C       ========================================================
C       Purpose : Compute the zeros of Hermite polynomial Ln(x)
C                 in the interval [-inf, inf], and the corresponding
C                 weighting coefficients for Gauss-Hermite
C                 integration
C       Input :   n    --- Order of the Hermite polynomial
C                 X(n) --- Zeros of the Hermite polynomial
C                 W(n) --- Corresponding weighting coefficients
C       ========================================================
C
C       Copyright: Professor Jianming Jin
C       Department of Electrical and Computer Engineering
C       University of Illinois at Urbana-Champaign
C       461 William L Everitt Laboratory
C       1406 West Green Street
C       Urbana, IL 61801-2991
C
C     Rewritten by Göran Broström to accommodate for the 
C     'sqrt(2*pi)' factor in the standard normal density.
C     Not any more! (2006-08-04)
C
C     Updated to more Fortran 77 than 66! (2006-08-04)

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER n, modified
      DIMENSION X(N), W(N)
C     LOGICAL modified

C     To silence the compiler:
      z = 0.0
      hf = 0.0
      hd = 0.0
      
      HN=1.0D0/N
      ZL=-1.1611D0+1.46D0*N**0.5
      DO NR = 1, N / 2
         IF (NR.EQ.1) then
            Z=ZL
         else
            Z=Z-HN*(N/2+1-NR)
         endif
         IT=0
C         do while(IT .LE. 40 .AND. DABS((Z-Z0)/Z) .GT. 1.0D-15)
C     Rather: repeat ... until... must check!
 10      IT=IT+1
            Z0=Z
            F0=1.0D0
            F1=2.0D0*Z
            DO K = 2, N
               HF=2.0D0*Z*F1-2.0D0*(K-1.0D0)*F0
               HD=2.0D0*K*F1
               F0=F1
               F1=HF
            enddo
            P=1.0D0
            DO I = 1, NR - 1
               P = P * (Z - X(I))
            enddo
            FD=HF/P
            Q=0.0D0
            DO  I=1,NR-1
               WP=1.0D0
               DO J=1,NR-1
                  IF (J .NE. I) WP=WP*(Z-X(J))
               enddo
               Q=Q+WP
            enddo
            GD=(HD-Q*FD)/P
            Z=Z-FD/GD
C         enddo
         IF (IT.LE.40.AND.DABS((Z-Z0)/Z).GT.1.0D-15) GO TO 10
         X(NR)=Z
         X(N+1-NR)=-Z
         R=1.0D0
         DO K=1,N
            R=2.0D0*R*K
         enddo
         W(NR)=3.544907701811D0*R/(HD*HD)
         W(N+1-NR)=W(NR)
      enddo
      IF ( N .NE. 2 * INT(N/2) ) THEN
         R1=1.0D0
         R2=1.0D0
         DO J = 1, N
            R1=2.0D0*R1*J
            IF (J.GE.(N+1)/2) R2=R2*J
         enddo
         W(N/2+1)=0.88622692545276D0*R1/(R2*R2)
         X(N/2+1)=0.0D0
      ENDIF

C **** Change back!
C      do i = 1, n
C        w(i) = w(i) /  1.77245385091
C        x(i) = x(i) * 1.41421356237
C      enddo

C *** New change! (2006-08-04) ***
C *** A correction: Make it depend on 'modified'.
      if (modified .eq. 1) then
         do i = 1, n
            w(i) = w(i) * exp(x(i)**2)
         enddo
      endif
      
      RETURN
      END
