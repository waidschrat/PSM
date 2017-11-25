SUBROUTINE LTI_KALMAN_FULLA_WITHINPUT(LL,INFO,A,B,C,D,TIME,Y,U,SIG,S,X0,&
   dimN,dimY,dimU,dimX,DoseN,DoseTime,DoseAmt,DoseState,&
   XF,XP,PF,PP,YP,R,KGAIN)


! LTI KALMAN FILTER
! SUBROUTINE designed to be called from R
!
! INPUT GIVEN FROM R
!
!	Matrix		A 		[dimX , dimX]
!           B 		[dimX , dimU]
!			      C 		[dimY,  dimX]
!			      D		  [dimY , dimU]
!
!	Data		T		[dimN]
!			    Y		[dimY , dimN]
!			    U		[dimU , dimN]
!
! Model		 SIG	  [dimX, dimX]
 !			   S      [dimY , dimY]
 !			   X0	    [dimX , 1]
 !
 ! Dimensions:
 !                      dimN , dimY, dimU , dimX

!COM        Kalman filter for the linear time-invariant (LTI) model
!COM
!COM           dX = [A*X + B*U]*dt + SIG*dW
!COM            Y =  C*X + D*U + e, e~N(0,S)
!COM
!COM        Doses:
!COM        Dosing is only permitted at observation times but several doses
!COM        can occur at the same observation time
!COM          DoseN       - integer length of DoseVector
!COM          DoseTime    - vector double   [DoseN]
!COM          DoseAmt     - vector double   [DoseN]
!COM          DoseState   - vector integer  [DoseN]
!COM
!COM   CALLS
!COM        SUBROUTINE DEXPM
!COM        SUBROUTINE LINEAR_MODEL
!COM        SUBROUTINE DTRM
!

IMPLICIT NONE
!
!-----------------------------------------------------DEFINITIONS-----!
!
INTEGER,PARAMETER                   :: IDEG=7           ! Used in the Matrix Exponential
DOUBLE PRECISION,PARAMETER          :: PI=3.141592653589793D0


INTEGER,INTENT(IN) :: dimX,dimY,dimU,dimN
INTEGER,INTENT(IN) :: DoseN
DOUBLE PRECISION,DIMENSION(dimN),INTENT(IN) :: TIME
DOUBLE PRECISION,DIMENSION(dimY,dimN),INTENT(IN) :: Y
DOUBLE PRECISION,DIMENSION(dimU,dimN),INTENT(IN) :: U
DOUBLE PRECISION,DIMENSION(dimX,dimX),INTENT(IN) :: A
DOUBLE PRECISION,DIMENSION(dimX,dimU),INTENT(IN) :: B
DOUBLE PRECISION,DIMENSION(dimY,dimX),INTENT(IN) :: C
DOUBLE PRECISION,DIMENSION(dimY,dimU),INTENT(IN) :: D
DOUBLE PRECISION,DIMENSION(dimX,1),INTENT(IN) :: X0
DOUBLE PRECISION,DIMENSION(dimX,dimX),INTENT(IN) :: SIG
DOUBLE PRECISION,DIMENSION(dimY,dimY),INTENT(IN) :: S
! DOSING
DOUBLE PRECISION,DIMENSION(DoseN), INTENT(IN) :: DoseTime,DoseAmt
INTEGER,DIMENSION(DoseN), INTENT(IN) :: DoseState

! INTENT OUT
INTEGER,INTENT(OUT) :: INFO
DOUBLE PRECISION,INTENT(OUT) :: LL
DOUBLE PRECISION,DIMENSION(dimX,dimN),INTENT(OUT) :: XF,XP

DOUBLE PRECISION,DIMENSION(dimY,dimN),INTENT(OUT) :: YP
DOUBLE PRECISION,DIMENSION(dimX,dimX, dimN),INTENT(OUT) :: PF,PP
DOUBLE PRECISION,DIMENSION(dimY,dimY,dimN),INTENT(OUT)  :: R
DOUBLE PRECISION,DIMENSION(dimX,dimY,dimN),INTENT(OUT)  :: KGAIN

! INTERNALS
DOUBLE PRECISION :: TAU,DETR,PS
INTEGER :: I,J,K,PIVX(dimX),PIVY(dimY) !,IALLOC
LOGICAL :: DEBUG
LOGICAL, DIMENSION(dimY) :: LOBSINDEX
INTEGER :: YCOUNT,TOTALMISSINGOBS,ID

DOUBLE PRECISION,DIMENSION(dimY,dimN) :: YPERR
DOUBLE PRECISION,DIMENSION(dimY,dimY,dimN) :: RINV

DOUBLE PRECISION,DIMENSION(dimX,dimX) :: AINV,PHIS,EYEDIMX,SIGSIGT,TMPXX,INTL
DOUBLE PRECISION,DIMENSION(dimX,dimY) :: TMPXY , TMP2XY
DOUBLE PRECISION,DIMENSION(dimY,dimX) :: TMPYX , TMP2YX
DOUBLE PRECISION,DIMENSION(dimY,dimY) :: TMPYY , EYEDIMY , E , TMP2YY , TMP3YY
DOUBLE PRECISION,DIMENSION(dimX,dimU) :: TMPXU
DOUBLE PRECISION,DIMENSION(dimY) :: TMPY , ERR
DOUBLE PRECISION,DIMENSION(dimX) :: TMPX


!-------------------------------------------------INITIALIZATIONS-----!
!
DEBUG = .FALSE.
! DEBUG = .TRUE.

IF(DEBUG) THEN
 !PRINT * , "Debugging Mode enabled - Done within the FORTRAN code"
 CALL INTPR("Debugging Mode enabled - Done within the FORTRAN code",53,0,0)
END IF

INFO = 0     !  INITIALIZE ERROR MESSAGE FLAG (= 0: NO ERROR)


PS          =  1.0D0    ! Initial covariance scaling
XF          =  0.0D0    ! state variables before filtering
PF          =  0.0D0    ! covariance matrix of XF
XP          =  0.0D0    ! one-step prediction state variables
PP          =  0.0D0    ! covariance matrix of XP

YP          =  0.0D0    ! output prediction
R           =  1.0D300  ! inverse covariance matrix of YP
RINV        =  0.0D0    ! inverse covariance matrix R

LL          =  0.0D0    ! minus log-likelihood value
KGAIN       =  0.0D0    ! kalman gain
DETR        =  0.0D0
AINV        =  0.0D0
TOTALMISSINGOBS = 0


TMPXX       =  0.0D0    ! help variable of dimension (NX,NX)
TMPXY       =  0.0D0    ! help variable of dimension (NX,NY)
TMP2XY      =  0.0D0    
TMPYX       =  0.0D0    ! help variable of dimension (NY,NX)
TMP2YX      =  0.0D0    
TMPYY       =  0.0D0    ! Tmp variable of dimension dimY, dimY
TMP2YY      =  0.0D0
TMP3YY      =  0.0D0
TMPXU       =  0.0D0
TMPY        =  0.0D0
TMPX        =  0.0D0
INTL        =  0.0D0

SIGSIGT     =  0.0D0    ! PIM*TRANSPOSE(PIM) of dimension (NX,NX)
PIVX        =  0        ! pivot elements of dimension dimX
PIVY        =  0        ! -- dimY

EYEDIMX     = 0.0D0
DO I=1,dimX
 EYEDIMX(I,I) = 1.0D0
END DO

EYEDIMY = 0.0D0
DO I=1,dimY
 EYEDIMY(I,I) = 1.0D0
END DO


!--------------------------------------PRE-FILTERING COMPUTATIONS-----!
!   Calculate SIG*SIGT - Writes result in SIGSIGT

!   F95  CALL DGEMM(SIG, SIG,SIGSIGT,'N','T',1.0D0,0.0D0)
    CALL DGEMM('N','T',dimX,dimX,dimX,1.0D0,SIG,dimX,SIG,dimX, &
       0.0D0,SIGSIGT,dimX)

!COM------------------------------------------------------------------!
!COM   INITIAL STATE COVARIANCE MATRIX PP(:,:,1) AT TIME 'K=1'        !
!COM------------------------------------------------------------------!

!COM   Eqn. (1.118) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.]!
!COM------------------------------------------------------------------!
TAU = TIME(2)-TIME(1)
CALL DEXPM(A,SIGSIGT,dimX,TAU,PHIS,INTL,IDEG,INFO)
IF (INFO /= 0) THEN
 INFO = 150
 RETURN
END IF
PP(:,:,1) = PS*INTL

! Set the initial value for the states
XP(:,1)  = X0(:,1)


IF(DEBUG) THEN
 !PRINT * , "Initial Statecovariance computed"
 CALL INTPR("Initial Statecovariance computed",32,0,0)
END IF

!COM------------------------------------------------------------------!
!COM  SINGULARITY OF A (NON-SINGULAR)
!COM------------------------------------------------------------------!

AINV = A
PIVX = 0
INFO = 0
IF(DEBUG) THEN
 CALL INTPR("A inverse",9,0,0)
END IF

CALL DGETRF(dimX,dimX,AINV,dimX,PIVX,INFO)

LU_FACTORIZATION_A: IF(INFO < 0) THEN
 INFO = 155
 RETURN
END IF LU_FACTORIZATION_A



!COM------------------------------------------------------------------!
!COM  CASE: A IS NON-SINGULAR  (SINGA == 0  ->  COMPUTE INVERSE OF A) !
!COM------------------------------------------------------------------!
!
CALL DGETRI(dimX,AINV,dimX,PIVX,TMPX,dimX,INFO)
IF(INFO /= 0) THEN
 INFO = 160
 RETURN
END IF

IF (DEBUG) THEN
 CALL INTPR("*** STARTING LOOP", 17,0,0)
END IF


!------------------------------------------------KALMAN FILTERING-----!
!
KALMANLOOP: DO K = 1, dimN
!
!-------------------------------------------------PREDICTION PART-----!
!

!IF(K > 3) THEN
!	DEBUG = .FALSE.
!END IF

IF(DEBUG) THEN
 CALL INTPR("-------------------",10,K,0)
 CALL INTPR("LOOP:", 5,K,1)
 CALL DBLEPR("TIME:",5,TIME(K),1)
END IF

!COM------------------------------------------------------------------!
!COM   OUTPUT PREDICTION YP(:,K) AT TIME 'K'
!COM   REFERENCE:
!COM   Eqn. (1.28) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.] !
!COM------------------------------------------------------------------!

IF(DEBUG) THEN
 CALL INTPR("***Prediction YP",16,0,0)
END IF

! Calculate Y = C*X + D*U
! Store C*X temporary in YP

! F95   CALL DGEMV(C,XP(:,K),YP(:,K),1.0D0,0.0D0,'N')
! CALL DGEMV('N',dimY,dimX,1.0D0,C,dimY,XP(:,K),dimX,0.0D0,YP(:,K),dimY)
CALL DGEMV('N',dimY,dimX,1.0D0,C,dimY,XP(:,K), 1 ,0.0D0,YP(:,K), 1 )

! F95  CALL DGEMV(D,U(:,K),YP(:,K),1.0D0,1.0D0,'N')
! CALL DGEMV('N',dimY,dimU,1.0D0,D,dimY,U(:,K),dimU,1.0D0,YP(:,K),dimY)
CALL DGEMV('N',dimY,dimU,1.0D0,D,dimY,U(:,K), 1 ,1.0D0,YP(:,K), 1)

IF(DEBUG) THEN
 CALL DBLEPR("YP :",3 ,YP(:,K),dimY)
END IF


!COM------------------------------------------------------------------!
!COM  MISSING OBSERVATIONS
!COM  Eqn. (1.29) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.]  !
!COM------------------------------------------------------------------!


IF(DEBUG) THEN
 CALL INTPR("*** Missing Observations",24,0,0)
 CALL DBLEPR(" Y : " ,2, Y(:,K),dimY)
END IF

DO I = 1,DIMY
 !LOBSINDEX(I) = (Y(I,K) < HUGE(1.0D0) )
 LOBSINDEX(I) = (Y(I,K) < 1.0D200 )
END DO

YCOUNT = COUNT(LOBSINDEX)

! UPDATE number of missing observations
TOTALMISSINGOBS = TOTALMISSINGOBS + dimY - YCOUNT

IF(DEBUG) THEN
 CALL INTPR("ObsINDEX:",9 , LOBSINDEX, dimY)
 CALL INTPR("YCOUNT :",8 , YCOUNT,1)
END IF

! Create E matrix
E = EYEDIMY
J = 1
DO I=1,dimY
 IF(LOBSINDEX(I) ) THEN
  E(J,:) = E(I,:)
  J = J+1
 END IF
END DO

DO I=J,dimY
 E(J,:) = 0.0D0
END DO

IF(DEBUG) THEN
 CALL DBLEPR("E-Matrix:",9 , E , dimY*dimY)
END IF

!COM------------------------------------------------------------------!
!COM  OUTPUT PREDICTION COVARIANCE MATRIX R(:,:,K)
!COM  Eqn. (1.29) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.]  !
!COM------------------------------------------------------------------!
!
! Calculate R= E * C * P *T(C)*T(E) + E*S*T(E)
! First TMP = C*P

IF(DEBUG) THEN
 CALL INTPR("*** R Calculation",17,0,0)
END IF

! First 	E*C
! F77 call dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
CALL DGEMM('N','N', YCOUNT, dimX, dimY, 1.0D0, E , dimY, C , dimY, &
   0.0D0, TMPYX, dimY)

IF(DEBUG) THEN
 CALL DBLEPR("C:",2,C,dimY*dimX)
 CALL DBLEPR("E*C:",4,TMPYX,dimY*dimX)
END IF

! Then (E*C)*P
CALL DGEMM('N','N', YCOUNT, dimX, dimX, 1.0D0, TMPYX , dimY , PP(:,:,K), dimX , &
   0.0D0, TMP2YX, dimY)
IF(DEBUG) THEN
 CALL DBLEPR("PP:",3,PP(:,:,K),dimX*dimX)
 CALL DBLEPR("(E*C)*PP:",9,TMP2YX,dimY*dimX)
END IF
      
! THEN (E*C*P)*T(C)
CALL DGEMM('N','T', YCOUNT, dimY, dimX, 1.0D0, TMP2YX , dimY , C , dimY , &
   0.0D0, TMPYY, dimY)
IF(DEBUG) THEN
 CALL DBLEPR("(E*C*P)*T(C):",13,TMPYY,dimY*dimY)
END IF
      
! THEN (E*C*P*T(C) ) * T(E)
CALL DGEMM('N','T', YCOUNT, YCOUNT, dimY, 1.0D0, TMPYY , dimY , E , dimY , &
   0.0D0, TMP2YY, dimY)
   
IF(DEBUG) THEN
 CALL DBLEPR("(E*C*P*T(C) )*T(E):",19,TMP2YY,dimY*dimY)
END IF

   
IF(DEBUG) THEN
 CALL DBLEPR("S:",2,S,dimY*dimY)
END IF
! THEN E*S
CALL DGEMM('N','N', YCOUNT, dimY, dimY, 1.0D0, E , dimY , S , dimY , &
   0.0D0, TMPYY, dimY)
IF(DEBUG) THEN
 CALL DBLEPR("E*S:",4,TMPYY,dimY*dimY)
END IF
   
! THEN (E*S) * T(E) 
CALL DGEMM('N','T', YCOUNT, YCOUNT, dimY, 1.0D0, TMPYY , dimY , E , dimY , &
   0.0D0, TMP3YY, dimY)

IF(DEBUG) THEN
 CALL DBLEPR("(E*S) * T(E):",13 , TMP3YY,dimY*dimY)
END IF
   


! Store in R
DO I=1,dimY
 DO J=1,dimY
  R(I,J,K) = TMP2YY(I,J) + TMP3YY(I,J)            
 END DO
END DO
   
IF(DEBUG) THEN
 CALL DBLEPR("R : ",3 , R(:,:,K),dimY*dimY)
END IF



! IF OBSERVATIONS ARE PRESENT NO UPDATING SHOULD BE PERFORMED
NO_OBS: IF(YCOUNT > 0) THEN


!COM------------------------------------------------------------------!
!COM  KALMAN GAIN KGAIN
!COM  Eqn. (1.31) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.]  !
!COM------------------------------------------------------------------!

IF(DEBUG) THEN
 CALL INTPR("*** Kalman gain",16,0,0)
END IF

! Calculate the Inverse R
PIVY = 0                ! PIVOT INDICIES FOR LU-FACTORIZATION
INFO = 0

TMPYY = R(:,:,K)

CALL DGETRF(YCOUNT,YCOUNT,TMPYY,dimY,PIVY,INFO)
IF(INFO < 0) THEN
 INFO = 155
 RETURN
ELSE IF(INFO > 0) THEN
 INFO = 156
 RETURN
END IF

CALL DGETRI(YCOUNT,TMPYY,dimY,PIVY,TMPY,dimY,INFO)
IF(INFO /= 0) THEN
 INFO = 160
 RETURN
END IF
!

! Remember  that only upper [YCOUNT ; YCOUNT] is valid
RINV(:,:,K) = TMPYY


IF(DEBUG) THEN
 CALL DBLEPR("RINV(:,:,K):",12, RINV(:,:,K),dimY*dimY)
END IF


! Calculate kalman gain
! K = P*T(C)*T(E)*Inv(R)

!   P*T(C)
CALL DGEMM('N','T', dimX, dimY, dimX, 1.0D0, PP(:,:,K), dimX, &
   C, dimY, 0.0D0, TMPXY, dimX)

! (P * T(C) ) * T(E)
CALL DGEMM('N','T',dimX,YCOUNT,dimY,1.0D0,TMPXY,dimX, &
   E, dimY, 0.0D0 , TMP2XY,dimX)
   
! (P *T(C)*T(E) ) * INV(R)
CALL DGEMM('N','N', dimX, YCOUNT, YCOUNT, 1.0D0, TMP2XY, dimX, &
   RINV(:,:,K), dimY, 0.0D0, KGAIN(:,:,K), dimX)

IF(DEBUG) THEN
 CALL DBLEPR("KGain : P*T(C)*T(E)*Inv(R) : ", 28 , KGAIN(:,:,K), dimX*dimX)
END IF


!------------------------------------------------------------------!
! LikeLihood contribution
!------------------------------------------------------------------!
! Add contribution to the loglikelihood from this observation

IF(DEBUG) THEN
 CALL INTPR("*** Likelihood Contribution",27,0,0)
END IF

! .5*(ln(det(R)) + e*Inv(R)*eT)
YPERR(:,K)  = Y(:,K) - YP(:,K)     ! INNOVATION AT TIME K


! Missing obs - No contribution
DO I=1,dimY
 IF(.NOT. LOBSINDEX(I)) THEN
  YPERR(I,K) = 0.0D0
 END IF
END DO


! Missing Observation -> fit ERR accordingly
! DGEMV(TRANSA, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
CALL DGEMV('N',dimY,dimY,1.0D0,E,dimY,YPERR(:,K),1,&
   0.0D0,ERR, 1)

IF(DEBUG) THEN
 CALL DBLEPR("ERR :"  ,5 , ERR,dimY)
END IF

CALL DGEMV('N',YCOUNT,YCOUNT,1.0D0,RINV(:,:,K),dimY,ERR,1,&
   0.0D0,TMPY,1)

! DETR <- Determinant of R
PIVY = 0
CALL DTRM( R(:,:,K),YCOUNT,DETR,PIVY)


! LL = LL  + LOG(DETR) + DOT(ERR,TMPY)
LL = LL  + 0.5D0*(LOG(DETR) + DOT_PRODUCT(ERR(1:YCOUNT),TMPY(1:YCOUNT)))

IF(DEBUG) THEN
 CALL DBLEPR("log(DETR) :",11,LOG(DETR),1)
 CALL DBLEPR("wE : ",4 ,  DOT_PRODUCT(ERR(1:YCOUNT),TMPY(1:YCOUNT)),1)
 CALL DBLEPR("dLL : ",5  , 0.5D0*(LOG(DETR) + DOT_PRODUCT(ERR(1:YCOUNT),TMPY(1:YCOUNT))) ,1)
END IF



!COM------------------------------------------------------------------!
!COM   UPDATED STATE VARIABLE VECTOR XF(:,K)                          !
!COM   Eqn. (1.32) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.] !
!COM------------------------------------------------------------------!

IF (DEBUG) THEN
 CALL INTPR("*** Update State variable XF=XP + KGAIN*E",41,0,0)
ENDIF

! XF = XP + KGAIN *E*YPERR

! DGEMV(TRANSA, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
CALL DGEMV('N',YCOUNT,dimY,1.0D0,E,dimY,YPERR(:,K),1, &
   0.0D0, TMPY, 1)


CALL DGEMV('N',dimX,YCOUNT,1.0D0,KGAIN(:,:,K),dimX,TMPY, 1,&
   0.0D0,XF(:,K), 1)


XF(:,K) = XP(:,K) + XF(:,K)

IF (DEBUG) THEN
 CALL DBLEPR("XF : ",4 , XF(:,K),dimX)
ENDIF


!
!COM------------------------------------------------------------------!
!COM   UPDATED STATE COVARIANCE MATRIX PF(:,:,K)                      !
!COM   Eqn. (1.33) [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.] !
!COM------------------------------------------------------------------!
!
! MATH FORMULAE
! PF = PP - K*R*T(K)


IF(DEBUG) THEN
 CALL INTPR("*** Updating State Covariance matrix PF= PP - K*R*T(K)",54,0,0)
END IF


CALL DGEMM('N','N',dimX, YCOUNT, YCOUNT, 1.0D0, KGAIN(:,:,K), dimX, &
   R(:,:,K), dimY, 0.0D0, TMPXY, dimX)

CALL DGEMM('N','T', dimX, dimX, YCOUNT, 1.0D0, TMPXY, dimX, &
   KGAIN(:,:,K), dimX, 0.0D0, TMPXX, dimX)


DO I=1,dimX
 DO J=1,dimX
  PF(I,J,K)=PP(I,J,K)-TMPXX(I,J)
 END DO
END DO


IF(DEBUG) THEN
 CALL DBLEPR("PF : ",4 , PF(:,:,K),dimX*dimX)
END IF

ELSE
 ! No OBSERVATION AT TIME[K]  -> No Update Only Pred
 XF(:,K) = XP(:,K)
 PF(:,:,K) = PP(:,:,K)

 IF(DEBUG) THEN
  CALL INTPR("No Observations -> Filtered = Prediction",40,0,0)
  CALL DBLEPR("XF : ",4 , XF(:,K),dimX)
  CALL DBLEPR("PF : ",4 , PF(:,:,K), dimX*dimX)
 END IF

END IF NO_OBS



!COM------------------------------------------------------------------!
!COM   DOSING
!COM------------------------------------------------------------------!

! DOSING Prediction and Updating finished - Dose if necesary
DO I=1,DoseN
 ! CHECK IF DOSING SHOULD OCCUR
 IF( DoseTime(I)==TIME(K) ) THEN

  IF(DEBUG) THEN
   CALL DBLEPR("Dosing at time " ,14, TIME(K) ,1)
   CALL INTPR("in CMT:",7,DoseState(I),1)
  ENDIF

  XF(DoseState(I),K) = XF(DoseState(I),K) + DoseAmt(I)

 END IF
END DO

!COM------------------------------------------------------------------!
!COM   NO MORE OBSERVATIONS
!COM------------------------------------------------------------------!

! IF THE OBSERVATIONS DATA IS FINISHED -> EXIT LOOP
IF(K==dimN) THEN
 EXIT
END IF

!COM------------------------------------------------------------------!
!COM   STATE PREDICTION COVARIANCE MATRIX PP(:,:,K+1)                 !
!COM   Eqn. (1.35) using Eqn. (1.45), (1.47), (1.48) and (1.49)       !
!COM   [CTSM 2.3 Math Guide, Dec. 2003, Kristensen, N.R.]             !
!COM------------------------------------------------------------------!

IF(DEBUG) THEN
 CALL INTPR("*** State Prediction Covariance matrix : PP",43,0,0)
END IF

IF(DEBUG) THEN
 CALL INTPR("Step1",5,0,0)
END IF


TAU = TIME(K+1)-TIME(K)

NO_TIMEDIFF: IF(TAU > 0.0D0) THEN

! DEXPM(A,SIGSIGT,DIMX,DT,H3T,H3TH2,IDEG,INFO)
INFO=0
CALL DEXPM(A,SIGSIGT,dimX,TAU,PHIS,INTL,IDEG,INFO)
IF (INFO /= 0) THEN
 INFO = 150
 RETURN
END IF


! F95 CALL DGEMM(PHIS, PF(:,:,K), PP(:,:,K+1), 'N', 'N', 1.0D0, 0.0D0)
CALL DGEMM('N','N', dimX, dimX, dimX, 1.0D0, PHIS, dimX, &
   PF(:,:,K), dimX, 0.0D0, TMPXX ,dimX)

! F95 CALL DGEMM(PP(:,:,K+1), PHIS, PP(:,:,K+1), 'N', 'T', 1.0D0, 0.0D0)
CALL DGEMM('N','T',dimX,dimX,dimX,1.0D0,TMPXX,dimX,PHIS,dimX,&
   0.0D0,PP(:,:,K+1),dimX)


DO I=1,dimX
 DO J=1,dimX
  PP(I,J,K+1) = PP(I,J,K+1) + INTL(I,J)
 ENDDO
ENDDO

IF(DEBUG) THEN
 CALL DBLEPR("PP : ",4 , PP(:,:,K+1),dimX*dimX)
END IF



!COM------------------------------------------------------------------!
!COM   STATE PREDICTION XP(:,K+1)                                     !
!COM------------------------------------------------------------------!
!COM    SPECIAL CASE NO. 3: NON-SINGULAR A, ZERO-ORDER HOLD ON INPUTS !
!COM    Eqn. (1.65) and (1.66) [CTSM 2.3 Math Guide, Dec. 2003,       !
!COM    Kristensen, N.R.]                                             !
!COM------------------------------------------------------------------!
!
! First Part of (1.65)

! XP = PHI*XF + AINV*(PHIS-EYE)*B*U

IF(DEBUG) THEN
 CALL INTPR("*** State Prediction  XP = PHI*XF + AINV*(PHIS-EYE)*B*U",55,0,0)
ENDIF

! F95 CALL DGEMV(PHIS, XF(:,K), XP(:,K+1), 1.0D0,0.0D0,'N')
CALL DGEMV('N',dimX,dimX,1.0D0,PHIS,dimX,XF(:,K), 1 ,&
   0.0D0,XP(:,K+1), 1 )


! F95 CALL DGEMM(AINV,PHIS-EYEDIMX,TMPXX,'N','N',1.0D0,0.0D0)
CALL DGEMM('N','N',dimX,dimX,dimX,1.0D0,AINV,dimX,PHIS-EYEDIMX,dimX,&
   0.0D0,TMPXX,dimX)


! F95 CALL DGEMM(TMPXX, B, TMPXU, 'N', 'N', 1.0D0, 0.0D0)
CALL DGEMM('N','N', dimX,dimU,dimX,1.0D0, TMPXX,dimX,B,dimX,&
   0.0D0,TMPXU,dimX)


! F95 CALL DGEMV(TMPXU,U(:,K),XP(:,K+1),1.0D0,1.0D0,'N')
CALL DGEMV('N',dimX,dimU,1.0D0,TMPXU,dimX,U(:,K), 1 ,&
   1.0D0,XP(:,K+1), 1)

   ELSE
      ! No Timedifference between observations - predicted state = filtered
     PP(:,:,K+1) = PF(:,:,K)
     XP(:,K+1) = XF(:,K)

   END IF NO_TIMEDIFF

END DO KALMANLOOP

IF(DEBUG) THEN
 CALL INTPR("**** Loop finished",18,0,0)
ENDIF

!---------------------------------------------------------------------!
! COMPUTE MINUS LOG-LIKELIHOOD VALUE: -LN(L(PHI,Y-Y0))                !
!---------------------------------------------------------------------!
!Update LL with dimY and dimN - NEGATIVE INDIVIDUAL LOG-LIKELIHOOD
LL = LL + 0.5D0*(LOG(2.0D0*PI)*(dimN*dimY-TOTALMISSINGOBS) )

! Set No Errors
INFO = 0

IF(DEBUG) THEN
 CALL DBLEPR("Inside Fortran LL: " ,19 ,LL,1)
END IF

MISSINGOBS: DO K = dimN,1,-1
DO I = 1,DIMY
 LOBSINDEX(I) = (Y(I,K) < 1.0D200 )
END DO
   
ID = COUNT(LOBSINDEX)**2
YCOUNT = COUNT(LOBSINDEX)

! NAs in YP
IF( YCOUNT.GT.0) THEN
 DO I=1,dimY
  IF(.NOT. LOBSINDEX(I)) THEN
   YP(I,K) = 1.0D300
  END IF
 END DO
ELSE
 YP(:,K) = 1.0D300
END IF
   
! Unfolding R
IF (YCOUNT.GT.0) THEN
 DO J = DIMY,1,-1
  DO I = DIMY,1,-1
   IF (LOBSINDEX(I).AND.LOBSINDEX(J)) THEN
    R(I,J,K) = R( MOD(ID-1,YCOUNT) +1,  (ID-1)/YCOUNT +1 ,K)
    ID = ID-1
   ELSE
    R(I,J,K) = 1.0D300
   ENDIF
  END DO
 END DO
ELSE
 R(:,:,K) = 1.0D300
END IF

      
! Unfolding of Kalman Gain	   
IF (YCOUNT.GT.0) THEN
 ID = YCOUNT
 DO J=dimY,1,-1
  IF( LOBSINDEX(J) ) THEN
   DO I=1,dimX
    KGAIN(I,J,K) = KGAIN(I,ID,K)
   END DO
   ID = ID-1
  ELSE
   DO I=1,dimX
    KGAIN(I,J,K) = 1.0D300
   END DO
  END IF
 END DO
ELSE
 KGAIN(:,:,K) = 1.0D300
END IF
   
   
END DO MISSINGOBS
!---------------------------------------END SUBROUTINE LTI_KALMAN-----!
END SUBROUTINE LTI_KALMAN_FULLA_WITHINPUT
