!---------------------------------------------------------------------!
!                                                                     !
!     V A R I O U S    M A T H E M A T I C A L    U T I L I T I E S   !
!                                                                     !
!                       FORTRAN 95 PROCEDURES                         !
!---------------------------------------------------------------------!
!
SUBROUTINE DTRM (A,N,D,INDX)
!---------------------------------------------------------------------!
!                                                                     !
! SUBROUTINE FOR EVALUATING THE DETERMINANT OF A MATRIX USING         !
! THE PARTIAL-PIVOTING GAUSSIAN ELIMINATION SCHEME.                   !
! COPYRIGHT (C) TAO PANG 2001.                                        !
!                                                                     !
! SOURCE: HTTP://WWW.PHYSICS.UNLV.EDU/~PANG/COMP3/CODE42.F90          !
!                                                                     !
!---------------------------------------------------------------------!
!                                                                     !
! A      - (N,N) MATRIX                                               !
! D      - DETERMINANT VALUE OF A                                     !
! INDX   - PIVOTING ORDER OF PARTIAL-PIVOTING GAUSS ELIMINATION METHOD!
!                                                                     !
! --------------------------------------------------------------------!

  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,MSGN
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  DOUBLE PRECISION, INTENT (OUT) :: D
  DOUBLE PRECISION, INTENT (INOUT), DIMENSION (N,N) :: A
!
  CALL ELGS(A,N,INDX)
!
  D = 1.0
  DO I = 1, N
    D = D*A(INDX(I),I)
  END DO
!
  MSGN = 1
  DO I = 1, N
    DO WHILE (I.NE.INDX(I))
          MSGN = -MSGN
          J = INDX(I)
          INDX(I) = INDX(J)
          INDX(J) = J
    END DO
  END DO
  D = MSGN*D
END SUBROUTINE DTRM

!
SUBROUTINE ELGS (A,N,INDX)
!
! SUBROUTINE TO PERFORM THE PARTIAL-PIVOTING GAUSSIAN ELIMINATION.
! A(N,N) IS THE ORIGINAL MATRIX IN THE INPUT AND TRANSFORMED MATRIX
! PLUS THE PIVOTING ELEMENT RATIOS BELOW THE DIAGONAL IN THE OUTPUT.
! INDX(N) RECORDS THE PIVOTING ORDER.  COPYRIGHT (C) TAO PANG 2001.
!
  IMPLICIT NONE
  INTEGER, INTENT (IN) :: N
  INTEGER :: I,J,K,ITMP
  INTEGER, INTENT (OUT), DIMENSION (N) :: INDX
  DOUBLE PRECISION :: C1,PI,PI1,PJ
  DOUBLE PRECISION, INTENT (INOUT), DIMENSION (N,N) :: A
  DOUBLE PRECISION, DIMENSION (N) :: C
!
! INITIALIZE THE INDEX
!
  DO I = 1, N
    INDX(I) = I
  END DO
!
! FIND THE RESCALING FACTORS, ONE FROM EACH ROW
!
  DO I = 1, N
    C1= 0.0
    DO J = 1, N
      C1 = DMAX1(C1,DABS(A(I,J)))
    END DO
    C(I) = C1
  END DO
!
! SEARCH THE PIVOTING (LARGEST) ELEMENT FROM EACH COLUMN
!
  DO J = 1, N-1
    PI1 = 0.0
    DO I = J, N
      PI = ABS(A(INDX(I),J))/C(INDX(I))
      IF (PI.GT.PI1) THEN
        PI1 = PI
        K   = I
      ENDIF
    END DO
!
! INTERCHANGE THE ROWS VIA INDX(N) TO RECORD PIVOTING ORDER
!
    ITMP    = INDX(J)
    INDX(J) = INDX(K)
    INDX(K) = ITMP
    DO I = J+1, N
      PJ  = A(INDX(I),J)/A(INDX(J),J)
!
! RECORD PIVOTING RATIOS BELOW THE DIAGONAL
!
      A(INDX(I),J) = PJ
!
! MODIFY OTHER ELEMENTS ACCORDINGLY
!
      DO K = J+1, N
        A(INDX(I),K) = A(INDX(I),K)-PJ*A(INDX(J),K)
      END DO
    END DO
  END DO
!
END SUBROUTINE ELGS
!---------------------------------------------------------------------!

! SOREN KLIM copied from 
! Andreas mods.f95 

SUBROUTINE DEXPM(A,SIGSIGT,DIMX,DT,H3T,H3TH2,IDEG,INFO)
!------------------------------------------------SUBROUTINE DEXPM-----
!
!COM------------------------------------------------------------------
!COM   MATRIX EXPONENTIAL AND MATRIX INTEGRAL                         
!COM------------------------------------------------------------------
!COM                                                                  
!COM   NAME                                                           
!COM        'DEXPM' - Subroutine                                      
!COM                                                                  
!COM   PURPOSE                                                        
!COM        Computes matrix exponential and matrix integral based on  
!COM        Eqn. (1.47), (1.48) and (1.49) [CTSM 2.3 Math Guide, Dec. 
!COM        2003, Kristensen, N.R.].                                  
!COM                                                                  
!COM   F95 INTERFACE                                                  
!COM        SUBROUTINE DEXPM(A,SIGSIGT,DT,H3T,H3TH2,IDEG,INFO)             
!COM                                                                  
!COM        USE SUNPERF                                               
!COM        USE MOPARAMS                                              
!COM                                                                  
!COM        INTEGER,INTENT(IN) :: IDEG                                
!COM        DOUBLE PRECISION,INTENT(IN) :: A(:,:),SIGSIGT(:,:),DT     
             
!COM        INTEGER,INTENT(OUT) :: INFO                               
!COM        DOUBLE PRECISION,INTENT(OUT) :: H3T(NX,NX),H3TH2(NX,NX)             
!COM                                                                  
!COM   REFERENCE                                                      
!COM        Eqn. (1.47), (1.48) and (1.49) [CTSM 2.3 Math Guide, Dec. 
!COM        2003, Kristensen, N.R.].                                  
!COM                                                                  
!COM   CALLS                                                          !
!COM        Subroutine DGPADM                                         !
!COM                                                                  !
!COM   ARGUMENTS                                                      !
!COM        A (input)                                                 !
!COM             Matrix A of the linear (LTI) model of dimensions     !
!COM             (NX,NX).                                             !
!COM                                                                  !
!COM        SIGSIGT (input)                                               !
!COM             SIGMA*TRANSPOSE(SIGMA) matrix of dimensions (NX,NX). !
!COM                                                                  !
!COM        DT (input)                                                !
!COM             Time interval T(k)-T(k-1). Scalar.                   !
!COM                                                                  !
!COM        H3T (output)                                              !
!COM             Matrix exponential EXP(A*DT) computed by means of a  !
!COM             Pade approximation. Matrix of dimensions (NX,NX).    !
!COM             Reference: Eqn. (1.47) and (1.48) [CTSM 2.3 Math     !
!COM             Guide,Dec. 2003, Kristensen, N.R.].                  !
!COM                                                                  !
!COM        H3TH2 (output)                                             !
!COM             Integral of EXP(A*DT)*SIGSIGT*TRANSPOSE(EXP(A*DT)) from  !
!COM             zero to DT. Matrix of dimensions (NX,NX).            !
!COM             Reference: Eqn. (1.47) and (1.49) [CTSM 2.3 Math     !
!COM             Guide,Dec. 2003, Kristensen, N.R.].                  !
!COM                                                                  !
!COM        IDEG (input)                                              !
!COM             Pade approximation order used in DGPADM (order of 6  !
!COM             is recommended by author of DGPADM).                 !
!COM                                                                  !
!COM        INFO (output)                                             !
!COM             = 0:  successful exit.                               !
!COM             > 0:  unsuccessful exit, see subroutine 'ERRORSTAT'. !
!COM                                                                  !
!COM------------------------------------------------------------------!
!
!-------------------------------------------------LINK TO MODULES-----!
!
! USE R TO COMPILE    USE SUNPERF
    IMPLICIT NONE
!
!-------------------------------------------VARIABLE DECLARATIONS-----!
    !
    INTEGER,INTENT(IN)  :: DIMX
    INTEGER,INTENT(IN)  :: IDEG
    DOUBLE PRECISION,INTENT(IN)  :: A(DIMX,DIMX),SIGSIGT(DIMX,DIMX),DT

    INTEGER,INTENT(OUT) :: INFO
    DOUBLE PRECISION,INTENT(OUT) :: H3T(DIMX,DIMX),H3TH2(DIMX,DIMX)
    
    ! INTERNALS
    INTEGER             :: I,IEXPH,NS,M,MM !,J
    DOUBLE PRECISION    :: G2(DIMX,DIMX) !,F3(DIMX,DIMX)
    
    INTEGER, ALLOCATABLE :: IPIV(:)
    DOUBLE PRECISION,ALLOCATABLE :: H(:,:) , WSP(:)
    
    ! Set Values
    !  INTEGER                     :: M=2*DIMX,MM=M*M
    M = 2*DIMX
    MM=M*M

    !PRINT * , "-- In MATUTIL95.f95 -> DEXPM"
    !PRINT * , "M : " , M
    !PRINT * , "MM : " , MM
    !PRINT * , "Lexp : " , 4*MM+IDEG+1

    ALLOCATE( IPIV(M) , H(M,M) , WSP(4*MM+IDEG+1) )
    
!
!-------------------------------------------------INITIALIZATIONS-----!
!
!
!   INITIALIZE ERROR MESSAGE FLAG (= 0: NO ERROR)
    INFO = 0
!
    H3TH2             = 0.0D0
    H(1:DIMX,1:DIMX)     = -A
    H(DIMX+1:M,1:DIMX)   = 0.0D0
    H(1:DIMX,DIMX+1:M)   = SIGSIGT
    H(DIMX+1:M,DIMX+1:M) = TRANSPOSE(A)
!
!COM------------------------------------------------------------------!
!COM  MATRIX EXPONENTIAL USING PADE APPROXIMATION                     !
!COM------------------------------------------------------------------!
!
    CALL DGPADM(IDEG,M,DT,H,M,WSP,4*MM+IDEG+1,IPIV,IEXPH,NS,INFO)
    IF (INFO /= 0) THEN
      INFO = 150
      RETURN
    END IF
!
!COM------------------------------------------------------------------!
!COM                                                                  !
!COM   H3T (in reference also named 'PHI')                            !
!COM                                                                  !
!COM   REFERENCE:                                                     !
!COM   Eqn. (1.47) and (1.48) [CTSM 2.3 Math Guide, Dec. 2003,        !
!COM   Kristensen, N.R.].                                             !
!COM------------------------------------------------------------------!
!COM   COPY I'TH COLUMN VECTORS IN LOWER RIGHT QUARDRANT OF EXPH INTO !
!COM   I'TH ROWS OF H3T, THUS RENDERING THE TRANSPOSE OF THE EXPH'S   !
!COM   LOWER RIGHT QUARDRANT INTO H3T.                                !
!COM------------------------------------------------------------------!
!
    DO I = DIMX+1, M
      H3T(I-DIMX,:) = WSP(IEXPH+(I-1)*M+DIMX:IEXPH+I*M-1)
    END DO
!
!COM------------------------------------------------------------------!
!COM                                                                  !
!COM   INTEGRAL[EXPA * SIGMA*TRANSPOSE(SIGMA) * TRANSPOSE(EXPA)]      !
!COM                                                                  !
!COM   REFERENCE:                                                     !
!COM   Eqn. (1.49) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.] !
!COM------------------------------------------------------------------!
!COM   COPY I'TH COLUMN VECTORS IN UPPER RIGHT QUARDRANT OF EXPH INTO !
!COM   I'TH ROWS OF H3T, THUS RENDERING THE TRANSPOSE OF THE EXPH'S   !
!COM   LOWER RIGHT QUARDRANT INTO H3T.                                !
!COM------------------------------------------------------------------!
!
    DO I = DIMX+1, M
      G2(:,I-DIMX) = WSP(IEXPH+(I-1)*M:IEXPH+I*M-1-DIMX)
    END DO
    CALL DGEMM('N','N',DIMX,DIMX,DIMX,1.0D0,H3T,DIMX,G2,DIMX,0.0D0,H3TH2,DIMX)
!  
!--------------------------------------------END SUBROUTINE DEXPM-----!
  END SUBROUTINE DEXPM
