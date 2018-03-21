module MatInv
	use types
	use global_var
	implicit none
contains

!SUBROUTINE SLVUNSYM(A,NRMAX,NCMAX,N,ITERM)
!   ___________________________________________________________
!   														   
!   Solver for BANDED UNSYMMETRIC system of algebraic equations
!   ___________________________________________________________
!                                                              
!	IMPLICIT REAL*8 (A-H,O-Z)
!   DOUBLE PRECISION A(NRMAX,NCMAX)
!      CERO=1.0D-15
!      PARE=CERO**2
!      NBND=2*ITERM
!      NBM=NBND-1
      
!     Begin elimination of the lower left
!      DO 80 I=1,N
!      IF (DABS(A(I,ITERM)).LT.CERO) GO TO 10    ! if1
!      GO TO 20
!   10 IF (DABS(A(I,ITERM)).LT.PARE) GO TO 110   ! if2
!   20 JLAST=MIN0(I+ITERM-1,N)
!      L=ITERM+1
!      DO 40 J=I,JLAST
!      L=L-1
!      IF (DABS(A(J,L)).LT.PARE) GO TO 40        ! if3
!      B=A(J,L)
!	   IF (B.LT.CERO) GO TO 40                   ! Xuan add for unstrict band matrix 08/24/2017
!      DO 30 K=L,NBND
!   30 A(J,K)=A(J,K)/B
!      IF (I.EQ.N) GO TO 90                      ! if4
!   40 CONTINUE
!      L=0
!      JFIRST=I+1
!      IF (JLAST.LE.I) GO TO 80                  ! if5
!      DO 70 J=JFIRST,JLAST
!      L=L+1
!	   IF (A(J,K-L).LT.CERO) GO TO 70            !Xuan add for unstrict band matrix 08/24/2017
!      IF (DABS(A(J,ITERM-L)).LT.PARE) GO TO 70  ! if6
!      DO 50 K=ITERM,NBM
!   50 A(J,K-L)=A(J-L,K)-A(J,K-L)
!      A(J,NBND)=A(J-L,NBND)-A(J,NBND)
!      IF (I.GE.N-ITERM+1) GO TO 70              ! if7
!      DO 60 K=1,L
!   60 A(J,NBND-K)=-A(J,NBND-K)
!   70 CONTINUE
!   80 CONTINUE
!   90 L=ITERM-1                                                                 
!      DO 100 I=2,N                                                              
!      DO 100 J=1,L                                                              
!      IF (N+1-I+J.GT.N) GO TO 100               ! if8
!      A(N+1-I,NBND)=A(N+1-I,NBND)-A(N+1-I+J,NBND)*A(N+1-I,ITERM+J)
!  100 CONTINUE
!      RETURN                                                                    
!  110 WRITE (6,140) I,A(I,ITERM)
!      STOP 

!  140 FORMAT (/,2X,'Computation stopped in SLVUNSYM:',I5,E12.4) 
!END SUBROUTINE SLVUNSYM

!======================================================================================================
SUBROUTINE SLVUNSYM(A,NRMAX,NCMAX,N,ITERM)
    implicit none                                                        
    DOUBLE PRECISION A(NRMAX,NCMAX)
	INTEGER NRMAX, NCMAX, N, ITERM

	DOUBLE PRECISION CERO, PARE
    INTEGER JFIRST, JLAST
	INTEGER L, NBND, NBM
	DOUBLE PRECISION B
    CERO=1.0D-15
    PARE=CERO**2
    NBND=2*ITERM
    NBM=NBND-1
      
!   Begin elimination of the lower left
    DO I=1,N
		IF (DABS(A(I,ITERM)).LT.PARE) THEN        ! IF2
			PRINT *, 'Computation stopped in SLVUNSYM:',I,A(I,ITERM)
		ENDIF
  	  	JLAST=MIN0(I+ITERM-1,N)
        L=ITERM+1
        DO J=I,JLAST
            L=L-1
			IF (DABS(A(J,L)).GT.PARE) THEN        ! IF3
		        B=A(J,L)                                                                  
		        DO K=L,NBND                                                            
		            A(J,K)=A(J,K)/B   
				END DO
				IF (I.EQ.N) GO TO 90
			ENDIF                                                      
		END DO
        L=0                                                             
        JFIRST=I+1
		IF (JLAST.GT.I) THEN					  ! IF5
        DO J=JFIRST,JLAST
            L=L+1
			IF (DABS(A(J,ITERM-L)).GE.CERO) THEN  ! IF6
		    DO K=ITERM,NBM
		        A(J,K-L)=A(J-L,K)-A(J,K-L)
			END DO
		    A(J,NBND)=A(J-L,NBND)-A(J,NBND)
			IF (I.LT.N-ITERM+1) THEN              ! IF7
			DO K=1,L
		   	    A(J,NBND-K)=-A(J,NBND-K)
			END DO
			ENDIF
			ENDIF
		END DO
		ENDIF
	END DO
  90 L=ITERM-1                                                                 
    DO I=2,N                                                              
    DO J=1,L                                                              
		IF (N+1-I+J.LE.N) THEN                    ! IF8
		  	A(N+1-I,NBND)=A(N+1-I,NBND)-A(N+1-I+J,NBND)*A(N+1-I,ITERM+J)
		ENDIF
	END DO
	END DO

	RETURN 
	STOP
END SUBROUTINE SLVUNSYM

end module MatInv

