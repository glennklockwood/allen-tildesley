********************************************************************************
** FICHE F.36.  MONTE CARLO SIMULATION OF HARD LINES IN 2D                    **
** This FORTRAN code is intended to illustrate points made in the text.       **
** To our knowledge it works correctly.  However it is the responsibility of  **
** the user to test it, if it is to be used in a research application.        **
********************************************************************************

C    *******************************************************************
C    ** TWO SEPARATE PARTS: FORTRAN AND BASIC VERSIONS.               **
C    *******************************************************************



C    *******************************************************************
C    ** FICHE F.36 - PART A                                           **
C    ** FORTRAN PROGRAM FOR MONTE CARLO SIMULATION OF HARD LINES      **
C    *******************************************************************



        PROGRAM HLINES

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** MONTE CARLO SIMULATION OF HARD LINES IN 2D.                   **
C    **                                                               **
C    ** THIS PROGRAM SETS UP AN INITIAL CONFIGURATION OF HARD LINES   **
C    ** AND CONDUCTS AN ENTIRE SWEEP IN DENSITY FROM LOW TO HIGH      **
C    ** OR HIGH TO LOW, PRINTING OUT THE ORDER PARAMETER REGULARLY.   **
C    ** THE LINE LENGTH IS TAKEN TO BE UNITY THROUGHOUT.              **
C    **                                                               **
C    ** PRINCIPAL VARIABLES:                                          **
C    **                                                               **
C    ** INTEGER N              NUMBER OF LINES                        **
C    ** REAL    RX(N),RY(N)    C-O-M POSITIONS OF LINES               **
C    ** REAL    EX(N),EY(N)    UNIT VECTOR ALONG LINE                 **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)

        INTEGER     MAXST
        PARAMETER ( MAXST = 15 )

        INTEGER     NEQUIL(MAXST), NPROD(MAXST)
        REAL        DENS(MAXST)

        REAL        MAXDEN
        PARAMETER ( MAXDEN = N / 4.0 )

        LOGICAL     LOSTRT
        INTEGER     STATE, NSTATE
        REAL        BOX, NEWBOX, MAXDIS, MAXROT, AVORD, RATDIS, RATROT

C    *******************************************************************

C    ** ENTER PARAMETERS **

        WRITE(*,'(1H1,'' **** PROGRAM HLINES ****                  '')')
        WRITE(*,'(//'' MONTE CARLO SIMULATION OF HARD LINES IN 2D  '')')
        WRITE(*,'('' HOW MANY STATE POINTS WOULD YOU LIKE ?        '')')
        READ (*,*) NSTATE
        WRITE(*,'('' LOW-DENSITY START CONFIGURATION (.T./.F.) ?   '')')
        READ (*,*) LOSTRT
        WRITE(*,'('' MAXIMUM ALLOWED DENSITY IS '',F10.6)') MAXDEN
        WRITE(*,'('' DENS   = DENSITY                              '')')
        WRITE(*,'('' NEQUIL = NUMBER OF SWEEPS FOR EQUILIBRATION   '')')
        WRITE(*,'('' NPROD  = NUMBER OF SWEEPS FOR PRODUCTION      '')')

        DO 10 STATE = 1, NSTATE

           WRITE(*,'('' DENS, NEQUIL, NPROD FOR STATE '',I5)') STATE
           READ (*,*) DENS(STATE), NEQUIL(STATE), NPROD(STATE)

           IF ( DENS(STATE) .GT. MAXDEN ) THEN

              WRITE(*,'('' SORRY - THAT ONE IS TOO BIG! '')')
              STOP

           ENDIF

10      CONTINUE

        WRITE(*,'('' NUMBER OF STATE POINTS       = '',I6)') NSTATE
        WRITE(*,'('' DESIRED PARAMETERS: '')')
        WRITE(*,'('' STATE  DENSITY  EQUILIBR   PRODUCT '')')

        DO 20 STATE = 1, NSTATE

           WRITE(*,'(1X,I5,F10.6,2I10)')
     +                 STATE, DENS(STATE), NEQUIL(STATE), NPROD(STATE)

20      CONTINUE

C    ** SET UP INITIAL CONFIGURATION **

        BOX = SQRT ( REAL ( N ) / DENS(1) )

        IF ( LOSTRT ) THEN

           WRITE(*,'('' SETTING UP LOW-DENSITY START '')')
           WRITE(*,'('' BOX LENGTH = '',F15.5)') BOX
           CALL LODEN ( BOX )

        ELSE

           WRITE(*,'('' SETTING UP HIGH-DENSITY START '')')
           WRITE(*,'('' BOX LENGTH = '',F15.5)') BOX
           CALL HIDEN ( BOX )

        ENDIF

        MAXDIS = 0.1
        MAXROT = 0.1

        WRITE(*,10001)

C    *******************************************************************
C    ** MAIN LOOP STARTS                                              **
C    *******************************************************************

        DO 1000 STATE = 1, NSTATE

           NEWBOX = SQRT( REAL(N) / DENS(STATE) )

           CALL EQUIL ( NEQUIL(STATE), BOX, NEWBOX,
     :                  MAXDIS, MAXROT, RATDIS, RATROT )

           DENS(STATE) = REAL(N) / BOX ** 2

           WRITE(*,10002)
     :        STATE, DENS(STATE), MAXDIS, RATDIS, MAXROT, RATROT

           CALL PRODUC ( NPROD(STATE), BOX, AVORD,
     :                   MAXDIS, MAXROT, RATDIS, RATROT )

           WRITE(*,10003)
     :        STATE, DENS(STATE), MAXDIS, RATDIS, MAXROT, RATROT, AVORD

1000    CONTINUE

C    *******************************************************************
C    ** MAIN LOOP ENDS                                                **
C    *******************************************************************

        STOP

10001   FORMAT (1X,'......STATE ...DENSITY ',
     :          '....MAXDIS RATDIS.... ....MAXROT RATROT.... ',
     :          'ORDER.....')
10002   FORMAT (1X,'EQUIL ',I5,1X,F10.5,1X,
     :           2(F10.5,1X,G10.3,1X))
10003   FORMAT (1X,'PRODU ',I5,1X,F10.5,1X,
     :           2(F10.5,1X,G10.3,1X),F10.5)

        END



        SUBROUTINE HIDEN ( BOX )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ROUTINE TO SET UP A HIGH DENSITY START CONFIGURATION          **
C    **                                                               **
C    ** LINE CENTRES ARE TRANSLATIONALLY DISORDERED, BUT ALL THE      **
C    ** LINES POINT IN THE SAME DIRECTION.  THIS DIRECTION IS CHOSEN  **
C    ** RANDOMLY.  THE CONFIGURATION CANNOT CONTAIN OVERLAPS, BUT WE  **
C    ** CHECK JUST IN CASE.                                           **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)
        REAL        BOX

        REAL        TWOPI
        PARAMETER ( TWOPI = 2.0 * 3.1415927 )

        INTEGER     I
        LOGICAL     OVRLAP
        REAL        TH, COSTH, SINTH, RANF, DUMMY, COSVAL, SINVAL

C    *******************************************************************

        TH = ( RANF ( DUMMY ) - 0.5 ) * TWOPI
        COSTH = COS ( TH )
        SINTH = SIN ( TH )

        DO 100 I = 1, N

           OVRLAP = .TRUE.

50         IF ( OVRLAP ) THEN

              RX(I) = ( RANF ( DUMMY ) - 0.5 ) * BOX
              RY(I) = ( RANF ( DUMMY ) - 0.5 ) * BOX
              EX(I) = COSTH
              EY(I) = SINTH
              OVRLAP = .FALSE.
              CALL DNCHEK ( BOX, I, OVRLAP )
              GOTO 50

           ENDIF

100     CONTINUE

        RETURN
        END



        SUBROUTINE LODEN ( BOX )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ROUTINE TO SET UP A LOW DENSITY START CONFIGURATION           **
C    **                                                               **
C    ** LINE CENTRES ARE TRANSLATIONALLY DISORDERED, AND ALL THE      **
C    ** LINE DIRECTIONS ARE SELECTED RANDOMLY.                        **
C    ** WE ENSURE NO OVERLAPS DURING CONSTRUCTION.                    **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)
        REAL        BOX

        REAL        TWOPI
        PARAMETER ( TWOPI = 2.0 * 3.1415927 )

        INTEGER     I
        LOGICAL     OVRLAP
        REAL        RANF, DUMMY, TH

C    *******************************************************************

        DO 100 I = 1, N

           OVRLAP = .TRUE.

50         IF ( OVRLAP ) THEN

              RX(I) = ( RANF ( DUMMY ) - 0.5 ) * BOX
              RY(I) = ( RANF ( DUMMY ) - 0.5 ) * BOX
              TH    = ( RANF ( DUMMY ) - 0.5 ) * TWOPI
              EX(I) = COS ( TH )
              EY(I) = SIN ( TH )
              OVRLAP = .FALSE.
              CALL DNCHEK ( BOX, I, OVRLAP )
              GOTO 50

           ENDIF

100     CONTINUE

        RETURN
        END



        SUBROUTINE EQUIL ( NSWEEP, BOX, NEWBOX,
     :                     MAXDIS, MAXROT, RATDIS, RATROT )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** CARRIES OUT NSWEEP SWEEPS OF EQUILIBRATION.                   **
C    **                                                               **
C    ** ATTEMPTS TO CHANGE BOX SIZE FROM BOX TO NEWBOX THROUGHOUT.    **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)

        INTEGER     NSWEEP
        REAL        BOX, NEWBOX, MAXDIS, MAXROT, RATDIS, RATROT

        LOGICAL     SCALED, OVRLAP
        INTEGER     SWEEP
        REAL        FACBOX, TRYBOX

C    *******************************************************************

        FACBOX = 1.0 / 2.0
        SCALED = .FALSE.

        DO 1000 SWEEP = 1, NSWEEP

C       ** CONDUCT USUAL MONTE CARLO MOVE SWEEP **

           CALL MOVE ( BOX, MAXDIS, MAXROT, RATDIS, RATROT )

C       ** NOW ATTEMPT TO SCALE BOX **

           IF ( .NOT. SCALED ) THEN

C          ** FIRST ATTEMPT COMPLETE SCALING **

              CALL SCALE ( BOX, NEWBOX )

              CALL CHECK ( NEWBOX, OVRLAP )

              IF ( OVRLAP ) THEN

C             ** IT FAILED: UNDO SCALING  **

                 SCALED = .FALSE.
                 CALL SCALE ( NEWBOX, BOX )

C             ** MAKE ONE ATTEMPT AT PARTIAL SCALING **

                 TRYBOX = BOX + FACBOX * ( NEWBOX - BOX )
                 CALL SCALE ( BOX, TRYBOX )
                 CALL CHECK ( TRYBOX, OVRLAP )

                 IF ( OVRLAP ) THEN

C                ** IT FAILED: UNDO SCALING AND   **
C                ** DECREASE FACBOX FOR NEXT TIME **

                    CALL SCALE ( TRYBOX, BOX )
                    FACBOX = FACBOX / 2.0

                 ELSE

C                ** IT WORKED: STICK WITH IT AND  **
C                ** INCREASE FACBOX FOR NEXT TIME **

                    BOX = TRYBOX
                    IF ( FACBOX .LT. 0.49 ) FACBOX = FACBOX * 2.0

                 ENDIF

              ELSE

C             ** IT WORKED FIRST TIME **

                 SCALED = .TRUE.
                 BOX = NEWBOX

              ENDIF

           ENDIF

1000    CONTINUE

        IF ( .NOT. SCALED ) THEN

           WRITE(*,'('' **** FAILED TO ACHIEVE NEW DENSITY ****'')')

        ENDIF

        NEWBOX = BOX

        RETURN
        END



        SUBROUTINE PRODUC ( NSWEEP, BOX, AVORD,
     :                      MAXDIS, MAXROT, RATDIS, RATROT )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** CARRIES OUT NSWEEP SWEEPS OF PRODUCTION.                      **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)

        INTEGER     NSWEEP
        REAL        BOX, AVORD, MAXDIS, MAXROT, RATDIS, RATROT

        INTEGER     SWEEP
        REAL        ORD

C    *******************************************************************

        AVORD = 0.0

        DO 1000 SWEEP = 1, NSWEEP

           CALL MOVE ( BOX, MAXDIS, MAXROT, RATDIS, RATROT )
           CALL ORDER ( ORD )
           AVORD = AVORD + ORD

1000    CONTINUE

        AVORD = AVORD / REAL ( NSWEEP )

        RETURN
        END



        SUBROUTINE MOVE ( BOX, MAXDIS, MAXROT, RATDIS, RATROT )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ATTEMPTS ONE TRANSLATIONAL OR ONE ROTATIONAL MOVE PER LINE    **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)
        REAL        MAXROT, MAXDIS, BOX, RATDIS, RATROT

        REAL        TWOPI
        PARAMETER ( TWOPI = 2.0 * 3.1415927 )

        LOGICAL     ROTMOV, OVRLAP
        INTEGER     I, TRYROT, TRYDIS, MAKROT, MAKDIS
        REAL        RANF, DUMMY, TH, COSTH, SINTH, DRX, DRY
        REAL        RXOLD, RYOLD, EXOLD, EYOLD
        REAL        MINDIS, MINROT

        PARAMETER ( MINDIS = 0.002, MINROT = 0.002 )

C    *******************************************************************

        TRYROT = 0
        TRYDIS = 0
        MAKROT = 0
        MAKDIS = 0

C    ** LOOP OVER ALL LINES **

        DO 1000 I = 1, N

           ROTMOV = ( RANF ( DUMMY ) .LT. 0.5 )

           IF ( ROTMOV ) THEN

C          ** ATTEMPT ROTATIONAL MOVE **

              TRYROT = TRYROT + 1
              TH     = ( RANF ( DUMMY ) - 0.5 ) * MAXROT
              COSTH  = COS ( TH )
              SINTH  = SIN ( TH )
              EXOLD  = EX(I)
              EYOLD  = EY(I)
              EX(I)  = EXOLD * COSTH - EYOLD * SINTH
              EY(I)  = EYOLD * COSTH + EXOLD * SINTH

              OVRLAP = .FALSE.
              CALL DNCHEK ( BOX, I, OVRLAP )
              CALL UPCHEK ( BOX, I, OVRLAP )

              IF ( OVRLAP ) THEN

                 EX(I) = EXOLD
                 EY(I) = EYOLD

              ELSE

                 MAKROT = MAKROT + 1

              ENDIF

           ELSE

C          ** ATTEMPT DISPLACEMENT **

              TRYDIS = TRYDIS + 1
              DRX    = ( RANF ( DUMMY ) - 0.5 ) * MAXDIS
              DRY    = ( RANF ( DUMMY ) - 0.5 ) * MAXDIS
              RXOLD  = RX(I)
              RYOLD  = RY(I)
              RX(I)  = RX(I) + DRX
              RY(I)  = RY(I) + DRY

              OVRLAP = .FALSE.
              CALL DNCHEK ( BOX, I, OVRLAP )
              CALL UPCHEK ( BOX, I, OVRLAP )

              IF ( OVRLAP ) THEN

                 RX(I) = RXOLD
                 RY(I) = RYOLD

              ELSE

                 MAKDIS = MAKDIS + 1

              ENDIF

           ENDIF

1000    CONTINUE

C    ** ADJUST MAXIMUM DISPLACEMENTS FOR NEXT TIME **

        IF ( TRYDIS .GT. 0 ) THEN

           RATDIS = REAL ( MAKDIS ) / REAL ( TRYDIS )

           IF ( RATDIS .GT. 0.5 ) THEN

              MAXDIS = MAXDIS + MINDIS

           ELSE

              MAXDIS = MAXDIS - MINDIS

           ENDIF

           IF ( MAXDIS .GT. BOX ) MAXDIS = BOX
           IF ( MAXDIS .LT. MINDIS ) MAXDIS = MINDIS

        ENDIF

        IF ( TRYROT .GT. 0 ) THEN

           RATROT = REAL(MAKROT) / REAL(TRYROT)

           IF ( RATROT .GT. 0.5 ) THEN

              MAXROT = MAXROT + MINROT

           ELSE

              MAXROT = MAXROT - MINROT

           ENDIF

           IF ( MAXROT .GT. TWOPI ) MAXROT = TWOPI
           IF ( MAXROT .LT. MINROT ) MAXROT = MINROT

        ENDIF

        RETURN
        END



        SUBROUTINE UPCHEK ( BOX, I, OVRLAP )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ROUTINE TO CHECK FOR OVERLAPS WITH LINES J>I                  **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)

        LOGICAL     OVRLAP
        INTEGER     I
        REAL        BOX

        REAL        TOL
        PARAMETER ( TOL = 1.0E-6 )

        INTEGER     J
        REAL        RXI, RYI, EXI, EYI
        REAL        RXIJ, RYIJ, RIJSQ, EXJ, EYJ, DET, DI, DJ
        REAL        BOXINV

C    *******************************************************************

        BOXINV = 1.0 / BOX

        RXI = RX(I)
        RYI = RY(I)
        EXI = EX(I)
        EYI = EY(I)

        J = I + 1

10      IF ( ( .NOT. OVRLAP ) .AND. ( J .LE. N ) ) THEN

           RXIJ = RXI - RX(J)
           RYIJ = RYI - RY(J)
           RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
           RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
           RIJSQ = RXIJ ** 2 + RYIJ ** 2

           IF ( RIJSQ .LT. 1.0 ) THEN

              EXJ = EX(J)
              EYJ = EY(J)

              DET = EXI * EYJ - EXJ * EYI

              IF ( ABS ( DET ) .GT. TOL ) THEN

                 DI = ( EXJ * RYIJ - EYJ * RXIJ ) / DET
                 DJ = ( EXI * RYIJ - EYI * RXIJ ) / DET

                 OVRLAP = ( ABS ( DI ) .LT. 0.5 ) .AND.
     :                    ( ABS ( DJ ) .LT. 0.5 )

              ENDIF

           ENDIF

           J = J + 1
           GOTO 10

        ENDIF

        RETURN
        END



        SUBROUTINE DNCHEK ( BOX, I, OVRLAP )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ROUTINE TO CHECK FOR OVERLAPS WITH LINES J<I                  **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)

        LOGICAL     OVRLAP
        INTEGER     I
        REAL        BOX

        REAL        TOL
        PARAMETER ( TOL = 1.0E-6 )

        INTEGER     J
        REAL        RXI, RYI, EXI, EYI
        REAL        RXIJ, RYIJ, RIJSQ, EXJ, EYJ, DET, DI, DJ
        REAL        BOXINV

C    *******************************************************************

        BOXINV = 1.0 / BOX

        RXI = RX(I)
        RYI = RY(I)
        EXI = EX(I)
        EYI = EY(I)

        J = I - 1

10      IF ( ( .NOT. OVRLAP ) .AND. ( J .GE. 1 ) ) THEN

           RXIJ = RXI - RX(J)
           RYIJ = RYI - RY(J)
           RXIJ = RXIJ - ANINT ( RXIJ * BOXINV ) * BOX
           RYIJ = RYIJ - ANINT ( RYIJ * BOXINV ) * BOX
           RIJSQ = RXIJ ** 2 + RYIJ ** 2

           IF ( RIJSQ .LT. 1.0 ) THEN

              EXJ = EX(J)
              EYJ = EY(J)

              DET = EXI * EYJ - EXJ * EYI

              IF ( ABS ( DET ) .GT. TOL ) THEN

                 DI = ( EXJ * RYIJ - EYJ * RXIJ ) / DET
                 DJ = ( EXI * RYIJ - EYI * RXIJ ) / DET

                 OVRLAP = ( ABS ( DI ) .LT. 0.5 ) .AND.
     :                    ( ABS ( DJ ) .LT. 0.5 )

              ENDIF

           ENDIF

           J = J - 1
           GOTO 10

        ENDIF

        RETURN
        END



        SUBROUTINE SCALE ( BOX, NEWBOX )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ROUTINE TO SCALE ALL COORDINATES FOR NEW BOX SIZE             **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)
        REAL        BOX, NEWBOX

        INTEGER     I
        REAL        FACTOR

C    *******************************************************************

        FACTOR = NEWBOX / BOX

        DO 1000 I = 1, N

           RX(I) = RX(I) * FACTOR
           RY(I) = RY(I) * FACTOR

1000    CONTINUE

        RETURN
        END



        SUBROUTINE CHECK ( BOX, OVRLAP )

C    *******************************************************************
C    ** ROUTINE TO CHECK ENTIRE CONFIGURATION FOR OVERLAPS            **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        BOX
        LOGICAL     OVRLAP

        INTEGER     I

C    *******************************************************************

        OVRLAP = .FALSE.
        I = 1

1000    IF ( ( .NOT. OVRLAP ) .AND. ( I .LT. N ) ) THEN

           CALL UPCHEK ( BOX, I, OVRLAP )
           I = I + 1
           GOTO 1000

        ENDIF

        RETURN
        END



        SUBROUTINE ORDER ( ORD )

        COMMON / BLOCK1 / RX, RY, EX, EY

C    *******************************************************************
C    ** ROUTINE TO COMPUTE ORDER PARAMETER                            **
C    **                                                               **
C    ** WE LOOK FOR THE MAXIMUM OF                                    **
C    ** SUM COS(2*TH(I) - 2*DI)                                       **
C    ** WHERE TH(I) IS THE ANGLE OF LINE I AND DI IS THE DIRECTOR     **
C    ** THIS OCCURS WHEN                                              **
C    ** (SIN(2*DI)/COS(2*DI)) = (SUM SIN(2*TH(I)))/(SUM COS(2*TH(I))) **
C    *******************************************************************

        INTEGER     N
        PARAMETER ( N = 25 )

        REAL        RX(N), RY(N), EX(N), EY(N)
        REAL        ORD

        INTEGER     I
        REAL        COS2TH, SIN2TH

C    *******************************************************************

        COS2TH = 0.0
        SIN2TH = 0.0

        DO 1000 I = 1, N

           COS2TH = COS2TH + ( 2.0 * EX(I) ** 2 - 1.0 )
           SIN2TH = SIN2TH + 2.0 * EX(I) * EY(I)

1000    CONTINUE

        ORD = SQRT ( SIN2TH ** 2 + COS2TH ** 2 ) / REAL ( N )

        RETURN
        END



C    *******************************************************************
C    ** FICHE F.36 - PART B                                           **
C    ** BASIC PROGRAM FOR MONTE CARLO SIMULATION OF HARD LINES        **
C    *******************************************************************



    1 MODE1
    2 CLS:PRINT''''TAB(12);"Program Hard Lines"
    3 PRINT'TAB(10);"Monte Carlo Simulation"
    4 PRINT'TAB(12);"of 2-D Hard Lines"
    5 PRINT''''TAB(15);"Fiche 36"
    6 INPUTTAB(10,31)"( press <RETURN> )" A$
   10DIM X(49),Y(49),T(49),CA(49),SA(49),XS1(49),YS1(49),XS2(49),YS2(49)
   15MODE4
   16@%=&20208
   20N=49
   30NR=7
   40INPUT "BOX LENGTH (400-700)";LB
   50LB2=LB/2
   60INPUT "LINE LENGTH ";L
   61INPUT "NUMBER OF CYCLES (1000)"; CYCLE
   62PRINT'''TAB(5);"BLEEPS FOR TRIAL OVERLAP"
   63 PRINT'TAB(5);"OUTPUTS ORDER PARAMETER"
   64 PRINTTAB(5);"EVERY FOUR CYCLES"
   68 INPUTTAB(5,31)"( press <RETURN> )" A$
   69LSQ=L*L
   70L2=L/2
   80DP=0
   90NA=0
  100DMAX=50
  110TMAX=PI/4
  120KSET=4*N
  130MAXKB=CYCLE*N
  134MODE1
  136COLOUR2
  137GCOL0,1
  140PROCBOX
  150PROCSETUP
  160I=0
  170REPEAT
  180KB=KB+1
  190I=I+1
  200IF I>N THEN I=1
  210DX=(2*RND(1)-1)*DMAX
  220DY=(2*RND(1)-1)*DMAX
  230DT=(2*RND(1)-1)*TMAX
  240XI=X(I)+DX
  250YI=Y(I)+DY
  260TI=T(I)+DT
  270IF XI>LB THEN XI=XI-LB ELSE IF XI<0 THEN XI=XI+LB
  280IF YI>LB THEN YI=YI-LB ELSE IF YI<0 THEN YI=YI+LB
  290IF TI>PI THEN TI=TI-PI ELSE IF TI<0 THEN TI=TI+PI
  300CI=COS(TI)
  310SI=SQR(1-CI^2)
  320PROCMOVE
  325K=0
  330K=K+1
  340IF K=I THEN 465
  350XDIFF=X(K)-XI
  360YDIFF=Y(K)-YI
  370XDIFF=XDIFF-LB*(XDIFF DIV LB2)
  380YDIFF=YDIFF-LB*(YDIFF DIV LB2)
  390RSQ=XDIFF^2+YDIFF^2
  400IF RSQ>LSQ THEN 465
  410S1=(CA(K)*SI-SA(K)*CI)*L
  420LM1=(YDIFF*CA(K)-XDIFF*SA(K))/S1
  430IF ABS(LM1)>0.5 THEN 465
  440LM2=(YDIFF*CI-XDIFF*SI)/S1
  450IF ABS(LM2)>0.5 THEN 465
  451SOUND 1,-10,53,5
  455PROCMOVEBACK
  460GOTO 540
  465REM
  470IF K<N THEN 330
  480X(I)=XI
  490Y(I)=YI
  500T(I)=TI
  510CA(I)=CI
  520SA(I)=SI
  521XS1(I)=XSI1
  522XS2(I)=XSI2
  523YS1(I)=YSI1
  524YS2(I)=YSI2
  530NA=NA+1
  540IF (KB MOD KSET)=0 THEN PROCDNEW
  550UNTIL KB=MAXKB
  560END
  565REM ***************************
  570DEF PROCDNEW
  575CBAR=0
  580RATIO=NA/KSET
  585TMAX=2*TMAX*RATIO
  590DMAX=2*DMAX*RATIO
  591@%=10
  592PRINT TAB(5,2);KB
  593@%=&20208
  596FORI=1 TO N
  597CBAR=CBAR+2*CA(I)^2-1
  598NEXT I
  599PRINT TAB(18,2);CBAR/N,DMAX
  600NA=0
  610ENDPROC
  615REM **************************
  620DEF PROCSETUP
  670FOR I=1 TO N
  700X(I)=LB*RND(1)
  710Y(I)=LB*RND(1)
  730NEXT I
  735TI=PI*RND(1)
  740FOR I =1 TO N
  750T(I)=TI
  760CA(I)=COS(TI)
  770SA(I)=SIN(TI)
  780XS1(I)=500-LB2+X(I)-L2*CA(I)
  790XS2(I)=XS1(I)+L*CA(I)
  800YS1(I)=500-LB2+Y(I)-L2*SA(I)
  810YS2(I)=YS1(I)+L*SA(I)
  820MOVE XS1(I),YS1(I)
  830PLOT 6,XS2(I),YS2(I)
  840NEXT
  850ENDPROC
  855REM ***************************
  860DEF PROCMOVE
  870XSI1=500-LB2+XI-L2*CI
  880XSI2=XSI1+L*CI
  890YSI1=500-LB2+YI-L2*SI
  900YSI2=YSI1+L*SI
  910MOVE XS1(I),YS1(I)
  920PLOT 6,XS2(I),YS2(I)
  930MOVE XSI1,YSI1
  940PLOT 6,XSI2,YSI2
  950ENDPROC
  955REM ***************************
  960DEF PROCMOVEBACK
  970MOVE XSI1,YSI1
  980PLOT 6,XSI2,YSI2
  990MOVE XS1(I),YS1(I)
 1000PLOT 6,XS2(I),YS2(I)
 1010ENDPROC
 1015REM **************************
 1020DEF PROCBOX
 1030MOVE 500-LB2,500-LB2
 1040DRAW 500+LB2,500-LB2
 1050DRAW 500+LB2,500+LB2
 1060DRAW 500-LB2,500+LB2
 1070DRAW 500-LB2,500-LB2
 1072PRINT TAB(1,1);"No. Configs.     O.P.      Max Disp."
 1080ENDPROC



